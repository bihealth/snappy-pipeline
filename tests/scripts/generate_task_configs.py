#!/usr/bin/env python3
"""Generate a task-based config with one task per workflow step and tool.

Each task starts from the step's `default_config_yaml()`. Three tables complete it:
- `TASK_CONFIG` fills config values and fixture paths per step and tool;
- `DEPENDENCY_TASKS` and `STEP_DEPENDENCY_TASKS` name the upstream task of each `depends_on` key,
  with per-tool rules in `wire_dependencies()`;
- `validate_and_autofill_step_config()` puts placeholders into required fields that remain unset.
"""

from __future__ import annotations

import argparse
import enum
import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, get_args, get_origin

import pydantic
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

yaml = ruamel_yaml.YAML()
yaml.preserve_quotes = True
yaml.default_flow_style = False
yaml.indent(sequence=4, offset=2)


#: Upstream task of each ``depends_on`` field.
DEPENDENCY_TASKS: dict[str, str] = {
    "alignments": "ngs_mapping_bwa",
    "annotated_variants": "variant_annotation_vep",
    "combined_variants": "combine_variants",
    "copy_number": "somatic_wgs_cnv_calling_cnvkit",
    "dbsnp": "dbsnp",
    "expression": "gene_expression_quantification_featurecounts",
    "features": "features",
    "fusions": "somatic_gene_fusion_calling_fusioncatcher",
    "germline_variants": "variant_filtration_bcftools",
    "hla_types": "hla_typing_optitype",
    "index": "reference_index_bwa",
    "phased_variants": "variant_phasing",
    "reads": "data_sets",
    "reference": "genome",
    "somatic_variants": "variant_calling_mutect2",
    "structural_variants": "sv_calling_targeted_delly2",
    "variants": "variant_calling_gatk4_hc_gvcf",
}

#: Exceptions to ``DEPENDENCY_TASKS``: (step, field) -> upstream task, or ``None`` to leave the
#: field unset. ``strandedness`` and ``panel_of_normals`` are wired per tool in build_all_tasks.
STEP_DEPENDENCY_TASKS: dict[tuple[str, str], str | None] = {
    # Expression export is off in the generated config, so it needs no RNA mapping task.
    ("cbioportal_export", "alignments"): None,
    ("cbioportal_export", "variants"): "variant_calling_mutect2",
    ("create_proteome", "variants"): "variant_annotation_vep",
    ("gene_expression_quantification", "alignments"): "ngs_mapping_star",
    ("reference_index", "reference"): "reference_download",
    ("homologous_recombination_deficiency", "copy_number"): (
        "somatic_targeted_seq_cnv_calling_sequenza"
    ),
    ("somatic_neoepitope_prediction", "germline_variants"): None,
    # Annotation of germline calls; revisit with the neoepitope backlog item in plans.md.
    ("somatic_neoepitope_prediction", "somatic_variants"): "variant_annotation_vep",
    ("somatic_targeted_seq_cnv_calling", "variants"): "variant_calling_mutect2",
    ("somatic_variant_signatures", "variants"): "variant_calling_mutect2",
    ("tumor_mutational_burden", "variants"): "variant_calling_mutect2",
    ("variant_export_external", "variants"): "external_vcf",
    ("variant_filtration", "variants"): "variant_annotation_vep",
    ("wgs_cnv_export_external", "variants"): "external_cnv",
    ("wgs_sv_export_external", "variants"): "external_sv",
    ("variant_phasing", "variants"): "variant_annotation_vep",
}


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("rt", encoding="utf-8") as f:
        data = yaml.load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping in {path}, got {type(data).__name__}")
    return data


def dump_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as f:
        yaml.dump(data, f)


def parse_default_step_config(workflow_cls: type, step_name: str) -> dict[str, Any]:
    text = (workflow_cls.default_config_yaml() or "").strip()
    if not text:
        return {}

    parsed = yaml.load(text)
    if not isinstance(parsed, dict):
        return {}

    step_config = parsed.get("step_config", {})
    if not isinstance(step_config, dict):
        return {}

    if step_name in step_config and isinstance(step_config[step_name], dict):
        return copy.deepcopy(step_config[step_name])

    return copy.deepcopy(step_config)


def _placeholder_for_annotation(annotation: Any) -> Any:
    origin = get_origin(annotation)
    args = get_args(annotation)

    if annotation in (str, Any):
        return "AUTO"
    if annotation is bool:
        return False
    if annotation is int:
        return 1
    if annotation is float:
        return 1.0
    if annotation is dict:
        return {}
    if annotation is list:
        return []
    if annotation is tuple:
        return []

    if origin in (list, tuple, set, frozenset):
        return []
    if origin is dict:
        return {}
    if origin is None and isinstance(annotation, type) and issubclass(annotation, enum.Enum):
        return next(iter(annotation)).value

    if origin is not None and args:
        # Handle Optional/Union by choosing the first non-None candidate.
        for candidate in args:
            if candidate is type(None):
                continue
            return _placeholder_for_annotation(candidate)

    if isinstance(annotation, type):
        # Nested Pydantic models and unknown classes can start as empty dicts.
        return {}

    return "AUTO"


def _set_nested(container: dict[str, Any], path: tuple[Any, ...], value: Any) -> None:
    current = container
    for key in path[:-1]:
        if not isinstance(key, str):
            return
        next_value = current.get(key)
        if not isinstance(next_value, dict):
            next_value = {}
            current[key] = next_value
        current = next_value

    last = path[-1]
    if isinstance(last, str):
        current[last] = value


def _unwrap_optional(annotation: Any) -> Any:
    origin = get_origin(annotation)
    args = get_args(annotation)
    if origin is not None and args:
        for candidate in args:
            if candidate is not type(None):
                return candidate
    return annotation


def _resolve_annotation_for_loc(config_model: type, loc: tuple[Any, ...]) -> Any:
    ann: Any = config_model
    for part in loc:
        if not isinstance(part, str):
            break
        ann = _unwrap_optional(ann)
        model_fields = getattr(ann, "model_fields", None)
        if not model_fields:
            break
        field_info = model_fields.get(part)
        if field_info is None:
            break
        ann = field_info.annotation
    return ann


def _remove_top_level_extra(config: dict[str, Any], path: tuple[Any, ...]) -> None:
    if path and isinstance(path[0], str):
        config.pop(path[0], None)


def _placeholder_for_error_type(err_type: str) -> Any:
    if err_type in {"string_type", "string_pattern_mismatch"}:
        return "AUTO"
    if err_type in {"int_type", "int_parsing"}:
        return 1
    if err_type in {"float_type", "float_parsing"}:
        return 1.0
    if err_type == "bool_type":
        return False
    if err_type == "list_type":
        return []
    if err_type == "dict_type":
        return {}
    if err_type in {"literal_error", "enum"}:
        return "AUTO"
    return "AUTO"


def validate_and_autofill_step_config(
    step_name: str,
    workflow_cls: type,
    step_config: dict[str, Any],
) -> tuple[dict[str, Any], list[str]]:
    notes: list[str] = []
    config_model = getattr(workflow_cls, "config_model_class", None)
    if not config_model:
        notes.append(f"{step_name}: no config_model_class, kept parsed defaults")
        return step_config, notes

    candidate = copy.deepcopy(step_config)
    max_rounds = 12
    for _ in range(max_rounds):
        try:
            validated = config_model.model_validate(candidate)
            return validated.model_dump(mode="json", exclude_none=True), notes
        except pydantic.ValidationError as e:
            changed = False
            for err in e.errors():
                err_type = err.get("type")
                loc = tuple(err.get("loc", ()))
                if not loc:
                    continue

                if err_type == "missing" and isinstance(loc[0], str):
                    resolved_ann = _resolve_annotation_for_loc(config_model, loc)
                    placeholder = _placeholder_for_annotation(resolved_ann)
                    _set_nested(candidate, loc, placeholder)
                    notes.append(
                        f"{step_name}: auto-filled missing field {'.'.join(map(str, loc))}"
                    )
                    changed = True
                elif err_type == "extra_forbidden":
                    _remove_top_level_extra(candidate, loc)
                    notes.append(f"{step_name}: removed extra field {'.'.join(map(str, loc))}")
                    changed = True
                elif err_type == "too_short" and isinstance(loc[0], str):
                    # Add a single-element placeholder list so min_length constraints are met.
                    resolved_ann = _resolve_annotation_for_loc(config_model, loc)
                    inner_args = get_args(_unwrap_optional(resolved_ann))
                    if inner_args:
                        inner_placeholder = _placeholder_for_annotation(inner_args[0])
                    else:
                        inner_placeholder = "AUTO"
                    _set_nested(candidate, loc, [inner_placeholder])
                    notes.append(
                        f"{step_name}: added single-item placeholder list for too_short field {'.'.join(map(str, loc))}"
                    )
                    changed = True
                elif isinstance(loc[0], str):
                    resolved_ann = _resolve_annotation_for_loc(config_model, loc)
                    placeholder = _placeholder_for_annotation(resolved_ann)
                    if placeholder in ({}, []) and err_type in {
                        "string_type",
                        "enum",
                        "literal_error",
                        "path_type",
                        "path_not_file",
                    }:
                        placeholder = _placeholder_for_error_type(err_type)
                    _set_nested(candidate, loc, placeholder)
                    notes.append(
                        f"{step_name}: coerced field {'.'.join(map(str, loc))} for error {err_type}"
                    )
                    changed = True

            if not changed:
                notes.append(f"{step_name}: unresolved validation errors, kept best-effort config")
                return candidate, notes
        except Exception as e:  # pragma: no cover - defensive for buggy validators
            notes.append(
                f"{step_name}: model validation raised {e.__class__.__name__}, kept best-effort config"
            )
            return candidate, notes

    notes.append(f"{step_name}: reached autofill iteration limit, kept best-effort config")
    return candidate, notes


def raw_data_dir(base_config: dict[str, Any], base_config_path: Path) -> str:
    data_sets = base_config.get("data_sets", {})
    if not isinstance(data_sets, dict) or not data_sets:
        return str(base_config_path.parent)

    first = next(iter(data_sets.values()))
    if not isinstance(first, dict):
        return str(base_config_path.parent)

    search_paths = first.get("search_paths", [])
    if isinstance(search_paths, list) and search_paths:
        return str((base_config_path.parent / str(search_paths[0])).resolve())

    return str(base_config_path.parent)


def dependency_fields(workflow_cls: type) -> list[str]:
    """Return the names of the ``depends_on`` fields of a step."""
    dep_field = workflow_cls.config_model_class.model_fields.get("depends_on")
    return [] if dep_field is None else list(dep_field.annotation.model_fields)


#: Test fixtures (index directories, BED files, gCNV models)
FIXTURES = Path(__file__).resolve().parents[1] / "snappy_pipeline" / "fixtures"

#: An existing file for config fields that must name one but whose content is never read
PLACEHOLDER_FILE = str(Path(__file__).resolve())


@dataclass(frozen=True)
class Fixture:
    """A path that ``fixture_paths()`` resolves for the base config."""

    name: str


@dataclass(frozen=True)
class Overwrite:
    """A value that replaces the default config's value instead of filling an unset one."""

    value: Any


def fixture_paths(base_config: dict[str, Any], base_config_path: Path) -> dict[str, Any]:
    """Return the fixture and reference paths that ``TASK_CONFIG`` refers to."""
    tasks = base_config.get("tasks", [])
    genome = next((task for task in tasks if task.get("name") == "genome"), None)
    reference = genome["config"]["files"]["fasta"] if genome else PLACEHOLDER_FILE

    def resolve(path: str) -> str:
        return path if path.startswith("/") else str((base_config_path.parent / path).resolve())

    def gcnv_models(library: str) -> list[dict[str, str]]:
        models = FIXTURES / "gcnv_models"
        return [
            {
                "library": library,
                "contig_ploidy": str(models / "ploidy_model"),
                "model_pattern": str(models / "call_model_*"),
            }
        ]

    return {
        "reference": resolve(reference),
        "star_index": str(FIXTURES / "star_index"),
        "mehari_db": str(FIXTURES / "mehari_db"),
        "cnvkit_targets": str(FIXTURES / "cnvkit_targets.bed"),
        "cnvkit_antitargets": str(FIXTURES / "cnvkit_antitargets.bed"),
        "gcnv_targeted": gcnv_models("default"),
        "gcnv_wgs": gcnv_models("wgs"),
        "placeholder_file": PLACEHOLDER_FILE,
        "raw_data_dir": [raw_data_dir(base_config, base_config_path)],
    }


DNA = "extraction_type == 'dna'"
RNA = "extraction_type == 'rna'"
REFERENCE = Fixture("reference")
PLACEHOLDER = Fixture("placeholder_file")
VARIANT_EXPORT_EXTERNAL = {
    "path_refseq_ser": PLACEHOLDER,
    "path_ensembl_ser": PLACEHOLDER,
    "path_db": PLACEHOLDER,
}

MELT_FILES = {"jar_file": PLACEHOLDER, "me_refs_path": PLACEHOLDER, "genes_file": PLACEHOLDER}

#: Config of the generated tasks, by (step, tool); (step, None) applies to every tool of a step.
#: Values fill keys that the default config leaves unset (missing, None, "", "AUTO", [] or
#: ["AUTO"]); ``Overwrite`` values replace whatever is there.
TASK_CONFIG: dict[tuple[str, str | None], dict[str, Any]] = {
    ("adapter_trimming", "bbduk"): {"bbduk": {"adapter_sequences": [PLACEHOLDER]}},
    ("external_data", None): {
        "produces": {"type": "raw"},
        "search_paths": Fixture("raw_data_dir"),
        "search_patterns": [
            {
                "left": r"(?P<readgroup>.+)\.R1\.fastq\.gz",
                "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz",
            }
        ],
    },
    ("gene_expression_quantification", None): {"library_selection": RNA},
    ("gene_expression_quantification", "dupradar"): {
        "dupradar": {"dupradar_path_annotation_gtf": PLACEHOLDER}
    },
    ("gene_expression_quantification", "rnaseqc"): {
        "rnaseqc": {"rnaseqc_path_annotation_gtf": PLACEHOLDER}
    },
    ("gene_expression_quantification", "strandedness"): {
        "strandedness": {"path_exon_bed": PLACEHOLDER}
    },
    ("gene_expression_quantification", "salmon"): {
        "salmon": {"path_index": Fixture("star_index"), "path_transcript_to_gene": PLACEHOLDER}
    },
    ("gene_expression_report", None): {"library_selection": RNA},
    ("hla_typing", "arcashla"): {"library_selection": RNA},
    ("ngs_data_qc", "picard"): {"picard": {"programs": ["CollectAlignmentSummaryMetrics"]}},
    ("ngs_mapping", None): {
        "target_coverage_report": {"enabled": False, "path_target_interval_list_mapping": []}
    },
    ("ngs_mapping", "bwa"): {"library_selection": DNA},
    ("ngs_mapping", "bwa_mem2"): {"library_selection": DNA},
    ("ngs_mapping", "minimap2"): {"library_selection": DNA},
    ("ngs_mapping", "mbcs"): {
        "library_selection": DNA,
        "mbcs": {"mapping_tool": "bwa"},
        "bwa": {},
        "bqsr": {"common_variants": REFERENCE},
    },
    ("ngs_mapping", "star"): {
        "library_selection": RNA,
        "strandedness": {"path_exon_bed": PLACEHOLDER, "strand": -1, "threshold": 0.85},
    },
    ("panel_of_normals", "cnvkit"): {"cnvkit": {"path_target": Overwrite("")}},
    ("panel_of_normals", "mutect2"): {"mutect2": {"germline_resource": REFERENCE}},
    ("panel_of_normals", "purecn"): {
        "purecn": {"path_bait_regions": REFERENCE, "path_normals_list": Overwrite("")}
    },
    ("reference_index", "star"): {"reference_molecule": Overwrite("rna")},
    ("repeat_expansion", None): {"repeat_catalog": REFERENCE, "repeat_annotation": REFERENCE},
    ("somatic_gene_fusion_calling", None): {"library_selection": RNA},
    ("somatic_gene_fusion_calling", "arriba"): {"arriba": {"path_index": Fixture("star_index")}},
    ("somatic_gene_fusion_calling", "defuse"): {
        "defuse": {"path_dataset_directory": Fixture("star_index")}
    },
    ("somatic_gene_fusion_calling", "hera"): {
        "hera": {"path_index": Fixture("star_index"), "path_genome": REFERENCE}
    },
    ("somatic_gene_fusion_calling", "jaffa"): {
        "jaffa": {"path_reference_files": Fixture("star_index")}
    },
    ("somatic_hla_loh_calling", None): {
        "path_hla_dat": PLACEHOLDER,
        "path_hla_fasta": PLACEHOLDER,
        "path_picard_dir": Fixture("star_index"),
    },
    ("somatic_gene_fusion_calling", "pizzly"): {
        "pizzly": {
            "kallisto_index": PLACEHOLDER,
            "transcripts_fasta": REFERENCE,
            "annotations_gtf": PLACEHOLDER,
        }
    },
    ("somatic_gene_fusion_calling", "star_fusion"): {
        "star_fusion": {"path_ctat_resource_lib": Fixture("star_index")}
    },
    ("cbioportal_export", None): {
        "copy_number_alteration": {"enabled": True},
        "path_gene_id_mappings": PLACEHOLDER,
    },
    ("somatic_msi_calling", None): {"loci_bed": REFERENCE},
    ("somatic_purity_ploidy_estimate", "ascat"): {"ascat": {"b_af_loci": PLACEHOLDER}},
    ("somatic_neoepitope_prediction", None): {
        "tool_hla_typing": {
            "dna": {"class_i": "optitype", "class_ii": None},
            "rna": {"class_i": None, "class_ii": None},
        }
    },
    ("somatic_targeted_seq_cnv_calling", "cnvkit"): {
        "cnvkit": {
            "path_target": Fixture("cnvkit_targets"),
            "path_antitarget": Fixture("cnvkit_antitargets"),
        }
    },
    ("somatic_targeted_seq_cnv_calling", "purecn"): {"purecn": {"path_container": PLACEHOLDER}},
    ("somatic_wgs_cnv_calling", "control_freec"): {
        "control_freec": {"path_chrlenfile": PLACEHOLDER, "path_mappability": PLACEHOLDER}
    },
    ("tumor_mutational_burden", None): {"target_regions": PLACEHOLDER},
    ("helper_gcnv_model_targeted", None): {"gcnv": {"path_uniquely_mapable_bed": PLACEHOLDER}},
    ("helper_gcnv_model_wgs", None): {"gcnv": {"path_uniquely_mapable_bed": PLACEHOLDER}},
    ("sv_calling_targeted", "melt"): {"melt": MELT_FILES},
    ("sv_calling_wgs", "melt"): {"melt": MELT_FILES},
    ("sv_calling_targeted", "gcnv"): {
        "gcnv": {"precomputed_model_paths": Fixture("gcnv_targeted")}
    },
    ("sv_calling_wgs", "gcnv"): {"gcnv": {"precomputed_model_paths": Fixture("gcnv_wgs")}},
    ("targeted_seq_mei_calling", "scramble"): {"scramble": {"blast_ref": REFERENCE}},
    ("variant_annotation", "mehari"): {
        "mehari": {"reference": REFERENCE, "transcripts": [PLACEHOLDER]}
    },
    ("variant_calling", None): {
        "baf_file_generation": {"enabled": False, "min_dp": 10},
        "bcftools_stats": {"enabled": False},
        "jannovar_stats": {"enabled": False, "path_ser": "AUTO"},
        "bcftools_roh": {
            "enabled": False,
            "path_af_file": "AUTO",
            "path_targets": None,
            "ignore_homref": False,
            "skip_indels": False,
            "rec_rate": 1e-8,
        },
    },
    ("varfish_export", None): {
        "path_mehari_db": Fixture("mehari_db"),
        "path_exon_bed": PLACEHOLDER,
    },
    ("variant_export_external", None): VARIANT_EXPORT_EXTERNAL,
    ("variant_phasing", None): {"gatk_read_backed_phasing": {"num_jobs": Overwrite(2)}},
    ("variant_filtration", "bcftools"): {"bcftools": {"exclude": "FILTER ~ 'low_depth'"}},
    ("variant_filtration", "regions"): {"regions": {"exclude": REFERENCE}},
    ("variant_filtration", "vembrane"): {"vembrane": {"expressions": {"some_filter": "True"}}},
    ("wgs_cnv_export_external", None): VARIANT_EXPORT_EXTERNAL,
    ("wgs_sv_export_external", None): VARIANT_EXPORT_EXTERNAL,
}


def _external_file(name: str, data_type: str, tags: list[str], key: str) -> dict[str, Any]:
    config = {"produces": {"type": data_type, "tags": tags}, "files": {key: REFERENCE}}
    return {"step": "external_data", "name": name, "config": config}


def _external_vcfs(name: str, tags: list[str]) -> dict[str, Any]:
    pattern = {"vcf": r".+\.vcf\.gz", "vcf_tbi": r".+\.vcf\.gz\.tbi"}
    config = {"produces": {"type": "variants", "tags": tags}, "search_patterns": [pattern]}
    config["search_paths"] = Fixture("raw_data_dir")
    return {"step": "external_data", "name": name, "config": config}


#: external_data tasks besides the step's own one: the reference data, and the inputs of the
#: external export steps. features and dbsnp name the reference FASTA, as no test reads them.
EXTRA_TASKS: list[dict[str, Any]] = [
    _external_file("genome", "raw", ["reference", "dna"], "fasta"),
    _external_file("features", "raw", ["features"], "gtf"),
    _external_file("dbsnp", "variants", ["dbsnp"], "vcf"),
    _external_vcfs("external_vcf", ["germline", "snv", "indel"]),
    _external_vcfs("external_cnv", ["germline", "cnv"]),
    _external_vcfs("external_sv", ["germline", "sv"]),
]


def _unset(value: Any) -> bool:
    return value is None or value in ("", "AUTO") or value in ([], ["AUTO"])


def _resolve(value: Any, fixtures: dict[str, Any]) -> Any:
    if isinstance(value, Fixture):
        return copy.deepcopy(fixtures[value.name])
    if isinstance(value, list):
        return [_resolve(item, fixtures) for item in value]
    if isinstance(value, dict):
        return {key: _resolve(item, fixtures) for key, item in value.items()}
    return value


def fill_config(config: dict[str, Any], fragment: dict[str, Any], fixtures: dict[str, Any]) -> None:
    """Fill the unset keys of ``config`` from ``fragment``; dicts are sections to recurse into."""
    for key, value in fragment.items():
        if isinstance(value, Overwrite):
            config[key] = _resolve(value.value, fixtures)
        elif isinstance(value, dict):
            if not isinstance(config.get(key), dict):
                config[key] = {}
            fill_config(config[key], value, fixtures)
        elif _unset(config.get(key)):
            config[key] = _resolve(value, fixtures)


def get_possible_tools(workflow_cls: type) -> list[str]:
    """Return the values of a step's ``tool`` enum; empty for steps without a tool."""
    tool_field = workflow_cls.config_model_class.model_fields.get("tool")
    if tool_field is None:
        return []
    annotation = tool_field.annotation
    while get_origin(annotation) not in (None, Literal):  # Annotated[...], X | None
        annotation = next(a for a in get_args(annotation) if a is not type(None))
    if get_origin(annotation) is Literal:
        return [str(value) for value in get_args(annotation)]
    if isinstance(annotation, type) and issubclass(annotation, enum.Enum):
        return [str(item.value) for item in annotation]
    return []


def wire_dependencies(step_name: str, tool: str | None) -> dict[str, str]:
    """Return the ``depends_on`` mapping of a generated task."""
    depends_on: dict[str, str] = {}
    for field in dependency_fields(WORKFLOW_REGISTRY[step_name]):
        upstream = STEP_DEPENDENCY_TASKS.get((step_name, field), DEPENDENCY_TASKS.get(field))
        if upstream:
            depends_on[field] = upstream

    # Expression quantifiers read the strandedness decision of the strandedness task.
    if step_name == "gene_expression_quantification" and tool in (
        "featurecounts",
        "dupradar",
        "duplication",
        "rnaseqc",
        "stats",
    ):
        depends_on["strandedness"] = "gene_expression_quantification_strandedness"

    # A mapper reads the index of its own tool (mbcs maps with bwa).
    if step_name == "ngs_mapping" and tool is not None:
        depends_on["index"] = f"reference_index_{'bwa' if tool == 'mbcs' else tool}"

    # The purecn panel of normals builds on the mutect2 one.
    if step_name == "panel_of_normals" and tool == "purecn":
        depends_on["panel_of_normals"] = f"{step_name}_mutect2"

    # cnvkit and purecn read the panel of normals built with the same tool.
    if step_name == "somatic_targeted_seq_cnv_calling" and tool in ("cnvkit", "purecn"):
        depends_on["panel_of_normals"] = f"panel_of_normals_{tool}"
    return depends_on


def build_all_tasks(base_config: dict[str, Any], base_config_path: Path) -> list[dict[str, Any]]:
    """Return one task per step and tool, with filled config and wired dependencies."""
    fixtures = fixture_paths(base_config, base_config_path)
    generation_notes: list[str] = []
    tasks: list[dict[str, Any]] = []
    for step_name, cls in sorted(WORKFLOW_REGISTRY.items()):
        for tool in get_possible_tools(cls) or [None]:
            config = parse_default_step_config(cls, step_name)
            if tool is not None:
                config["tool"] = tool
                # The selected tool's section must be present, even if empty.
                if tool in cls.config_model_class.model_fields:
                    config.setdefault(tool, {})
            fill_config(config, TASK_CONFIG.get((step_name, None), {}), fixtures)
            if tool is not None:
                fill_config(config, TASK_CONFIG.get((step_name, tool), {}), fixtures)
            if depends_on := wire_dependencies(step_name, tool):
                config["depends_on"] = depends_on

            config, notes = validate_and_autofill_step_config(step_name, cls, config)
            generation_notes.extend(notes)
            if depends_on:
                # Validation dumps every depends_on field; keep only the wired ones.
                config["depends_on"] = depends_on
            name = f"{step_name}_{tool}" if tool else step_name
            tasks.append({"step": step_name, "name": name, "config": config})

    tasks += [_resolve(extra, fixtures) for extra in EXTRA_TASKS]

    if generation_notes:
        print("Generation notes:")
        for note in generation_notes:
            print(f"- {note}")

    return tasks


def normalize_data_set_paths(config: dict[str, Any], base_config_path: Path) -> None:
    base_dir = base_config_path.parent
    data_sets = config.get("data_sets", {})
    if not isinstance(data_sets, dict):
        return

    for _, data_set in data_sets.items():
        if not isinstance(data_set, dict):
            continue

        sheet_file = data_set.get("file")
        if isinstance(sheet_file, str) and sheet_file:
            path = Path(sheet_file)
            if not path.is_absolute():
                data_set["file"] = str((base_dir / path).resolve())

        search_paths = data_set.get("search_paths", [])
        if isinstance(search_paths, list):
            normalized = []
            for p in search_paths:
                if isinstance(p, str) and p:
                    path = Path(p)
                    normalized.append(
                        str((base_dir / path).resolve()) if not path.is_absolute() else p
                    )
                else:
                    normalized.append(p)
            data_set["search_paths"] = normalized


def build_config(
    base_config: dict[str, Any], tasks: list[dict[str, Any]], base_config_path: Path
) -> dict[str, Any]:
    config = {"tasks": tasks, "data_sets": copy.deepcopy(base_config.get("data_sets", {}))}
    normalize_data_set_paths(config, base_config_path)
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=Path("tests/snappy_pipeline/fixtures/base_config.yaml"),
        help="Path to an existing config.yaml used as source for the genome task and data_sets.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("scratch/generated-configs"),
        help="Directory to write generated config folders to.",
    )
    args = parser.parse_args()

    base_config_path = args.base_config.resolve()
    base_config = load_yaml(base_config_path)

    all_tasks = build_all_tasks(base_config, base_config_path)

    all_config_dir = args.out_dir.resolve() / "all-workflows"
    all_config = build_config(base_config, all_tasks, base_config_path)
    dump_yaml(all_config_dir / "config.yaml", all_config)

    print(f"Wrote {all_config_dir / 'config.yaml'} with {len(all_tasks)} tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
