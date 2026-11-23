#!/usr/bin/env python3
"""Generate a task-based config with one task per workflow step and tool.

Each task starts from the step's `default_config_yaml()`. Three tables complete it:
- `TASK_CONFIG` fills config values and fixture paths per step and tool;
- `DEPENDENCY_TASKS` and `STEP_DEPENDENCY_TASKS` name the upstream task of each `depends_on` key,
  with per-tool rules in `wire_dependencies()`.

`validate_step_config()` drops the `depends_on` keys that no reader of the task reads;
any other validation error fails the generation.
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

from snappy_pipeline.models import DependsOnError
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

yaml = ruamel_yaml.YAML()
yaml.preserve_quotes = True
yaml.default_flow_style = False
yaml.indent(sequence=4, offset=2)


#: Upstream task of each ``depends_on`` field.
DEPENDENCY_TASKS: dict[str, str] = {
    "alignments": "ngs_mapping_bwa",
    "qc": "ngs_data_qc_alfred",
    "annotated_variants": "variant_annotation_vep",
    "candidates": "variant_calling_freebayes",
    "combined_variants": "combine_variants",
    "copy_number": "somatic_wgs_cnv_calling_cnvkit",
    "common_variants": "population_variants",
    "dbsnp": "dbsnp",
    "excluded_regions": "regions",
    "gene_regions": "tagged_regions",
    "germline_resource": "population_variants",
    "exons": "exons",
    "expression": "gene_expression_quantification_featurecounts",
    "features": "features",
    "fusions": "somatic_gene_fusion_calling_arriba",
    "germline_variants": "variant_filtration_bcftools",
    "hla_types": "hla_typing_optitype",
    "index": "reference_index_bwa",
    "kits": "kits",
    "mappable_regions": "mappable_regions",
    "microsatellites": "regions",
    "par_regions": "tagged_regions",
    "phased_variants": "variant_phasing_gatk_read_backed_phasing",
    "reads": "data_sets",
    "reference": "genome",
    "snp_loci": "tagged_regions",
    "somatic_variants": "variant_calling_mutect2",
    "structural_variants": "sv_calling_targeted_delly2",
    "variants": "variant_calling_gatk4_hc_gvcf",
}

#: Steps with one task per value of a config field instead of per tool: step -> field.
#: The tasks are named ``<step>_<value>``, and ``TASK_CONFIG`` keys them by the value.
VARIANT_FIELDS: dict[str, str] = {"gcnv_model": "intervals"}

#: Exceptions to ``DEPENDENCY_TASKS``: (step, field) -> upstream task, or ``None`` to leave the
#: field unset. ``strandedness`` and ``panel_of_normals`` are wired per tool in build_all_tasks.
STEP_DEPENDENCY_TASKS: dict[tuple[str, str], str | list[str] | None] = {
    ("cbioportal_export", "variants"): "variant_annotation_pvactools",
    # The small variants and the structural variants of each pedigree
    ("concat_variants", "variants"): [
        "variant_calling_varlociraptor_germline",
        "variant_calling_varlociraptor_sv",
    ],
    ("create_proteome", "variants"): "variant_annotation_vep",
    ("gene_expression_quantification", "alignments"): "ngs_mapping_star",
    ("gene_expression_quantification", "index"): "salmon_index",
    ("somatic_gene_fusion_calling", "index"): "star_index",
    ("reference_index", "reference"): "genome",
    ("homologous_recombination_deficiency", "copy_number"): (
        "somatic_targeted_seq_cnv_calling_sequenza"
    ),
    # Neoepitope prediction with every input: VEP with the pVACtools plugins, HLA types of DNA
    # and RNA, and the RNA alignments with their strandedness
    ("somatic_neoepitope_prediction", "variants"): "variant_annotation_pvactools",
    ("somatic_neoepitope_prediction", "hla_types"): ["hla_typing_optitype", "hla_typing_arcashla"],
    ("somatic_neoepitope_prediction", "alignments"): "ngs_mapping_star",
    ("somatic_neoepitope_prediction", "phased_variants"): "variant_phasing_pvactools",
    ("somatic_neoepitope_prediction", "expression"): "gene_expression_quantification_salmon",
    ("somatic_neoepitope_prediction", "fusions"): "somatic_gene_fusion_calling_arriba",
    ("somatic_neoepitope_prediction", "proteome"): "proteome",
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


def parse_default_step_config(workflow_cls: type) -> dict[str, Any]:
    config = yaml.load(workflow_cls.default_config_yaml())["tasks"][0]["config"]
    return copy.deepcopy(config) if isinstance(config, dict) else {}


def validate_step_config(name: str, workflow_cls: type, config: dict[str, Any]) -> dict[str, Any]:
    """Return the validated config of the generated task ``name``.

    Drops the ``depends_on`` fields that no reader of the task reads; any other error fails.
    """
    while True:
        try:
            validated = workflow_cls.config_model_class.model_validate(config)
        except pydantic.ValidationError as e:
            unread = [
                field
                for err in e.errors()
                if isinstance(error := err.get("ctx", {}).get("error"), DependsOnError)
                for field in error.unread
            ]
            if not unread:
                e.add_note(f"in the generated task {name}")
                raise
            for field in unread:
                config["depends_on"].pop(field)
        else:
            return validated.model_dump(mode="json", exclude_none=True)


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


#: Test fixtures (index directories, BED files)
FIXTURES = Path(__file__).resolve().parents[1] / "snappy_pipeline" / "fixtures"

#: An existing file for config fields that must name one but whose content is never read
PLACEHOLDER_FILE = str(Path(__file__).resolve())


@dataclass(frozen=True)
class Fixture:
    """A path that ``fixture_paths()`` resolves for the base config, optionally a path below it."""

    name: str
    subpath: str = ""


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

    return {
        "reference": resolve(reference),
        "star_index": str(FIXTURES / "star_index"),
        "cnvkit_targets": str(FIXTURES / "cnvkit_targets.bed"),
        "mappable": str(FIXTURES / "mappable.bed"),
        "mehari_db": str(FIXTURES / "mehari_db"),
        "placeholder_file": PLACEHOLDER_FILE,
        "vcf": str(FIXTURES.parents[2] / ".tests/test-workflow/resources/refs/dbsnp.vcf.gz"),
        "raw_data_dir": [raw_data_dir(base_config, base_config_path)],
    }


DNA = "extraction_type == 'dna'"
RNA = "extraction_type == 'rna'"
REFERENCE = Fixture("reference")
PLACEHOLDER = Fixture("placeholder_file")
VCF = Fixture("vcf")
MAPPABLE = Fixture("mappable")
VARIANT_EXPORT_EXTERNAL = {
    "path_refseq_ser": PLACEHOLDER,
    "path_ensembl_ser": PLACEHOLDER,
    "path_db": PLACEHOLDER,
}

MELT_FILES = {"jar_file": PLACEHOLDER, "me_refs_path": PLACEHOLDER}

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
    ("gene_expression_quantification", "dupradar"): {"dupradar": {}},
    ("gene_expression_quantification", "rnaseqc"): {"rnaseqc": {}},
    ("gene_expression_quantification", "salmon"): {"salmon": {}},
    ("ngs_data_qc", "picard"): {
        "picard": {"programs": ["CollectAlignmentSummaryMetrics", "CollectHsMetrics"]}
    },
    ("ngs_mapping", "bwa"): {"library_selection": DNA, "bam_collect_doc": {"enabled": True}},
    ("ngs_mapping", "bwa_mem2"): {"library_selection": DNA},
    ("ngs_mapping", "bwa_mem3"): {"library_selection": DNA},
    ("ngs_mapping", "minimap2"): {"library_selection": DNA},
    ("ngs_mapping", "star"): {
        "library_selection": RNA,
        "strandedness": {"strand": -1, "threshold": 0.85},
    },
    ("panel_of_normals", "purecn"): {"purecn": {"path_normals_list": Overwrite("")}},
    ("reference_download", None): {"source": "ensembl"},
    ("somatic_gene_fusion_calling", "arriba"): {"arriba": {}},
    ("somatic_hla_loh_calling", None): {
        "path_hla_dat": PLACEHOLDER,
        "path_hla_fasta": PLACEHOLDER,
        "path_picard_dir": Fixture("star_index"),
    },
    ("somatic_gene_fusion_calling", "star_fusion"): {"star_fusion": {}},
    ("cbioportal_export", None): {
        "copy_number_alteration": {"enabled": True},
        "path_gene_id_mappings": PLACEHOLDER,
        "vcf2maf": {"Center": "AUTO", "ncbi_build": "GRCh37"},
        "study": {
            "type_of_cancer": "AUTO",
            "cancer_study_id": "AUTO",
            "study_description": "AUTO",
            "study_name": "AUTO",
            "study_name_short": "AUTO",
        },
    },
    ("somatic_neoepitope_prediction", None): {
        "pvactools": {
            "pvacseq": {
                "pileup": {"enabled": Overwrite(True)},
                "expression": {"enabled": Overwrite(True)},
                "phasing": {"enabled": Overwrite(True)},
                "net_chop": {"enabled": Overwrite(True), "path": PLACEHOLDER},
                "netmhc_stab": {"enabled": Overwrite(True)},
            },
            "pvacfuse": {"enabled": Overwrite(True)},
            "pvacsplice": {"enabled": Overwrite(True)},
        }
    },
    ("somatic_wgs_cnv_calling", "control_freec"): {"control_freec": {"convert": {}}},
    ("sv_calling_targeted", "melt"): {"melt": {**MELT_FILES, "me_refs_infix": "AUTO"}},
    ("sv_calling_wgs", "melt"): {"melt": MELT_FILES},
    ("targeted_seq_mei_calling", "scramble"): {"scramble": {"blast_ref": REFERENCE}},
    ("variant_calling", None): {
        "baf_file_generation": {"enabled": False, "min_dp": 10},
        "bcftools_stats": {"enabled": False},
        "bcftools_roh": {
            "enabled": False,
            "path_af_file": "AUTO",
            "ignore_homref": False,
            "skip_indels": False,
            "rec_rate": 1e-8,
        },
    },
    # The statistics and their QC summary, once: one job per task, one file per pedigree member
    ("variant_calling", "bcftools_call"): {"bcftools_call": {"assembly": "GRCh37"}},
    ("variant_calling", "gatk4_hc_joint"): {"bcftools_stats": {"enabled": Overwrite(True)}},
    ("variant_calling", "freebayes"): {"freebayes": {"num_jobs": 2}},
    ("sv_calling_wgs", "popdel"): {"popdel": {"num_jobs": 2}},
    ("variant_calling", "mutect2"): {"mutect2": {"contamination": {}}},
    ("variant_calling", "varlociraptor"): {"varlociraptor": {"scenario": "somatic", "num_jobs": 2}},
    ("variant_export_external", None): VARIANT_EXPORT_EXTERNAL,
    ("variant_phasing", "gatk_read_backed_phasing"): {
        "gatk_read_backed_phasing": {"num_jobs": Overwrite(2)}
    },
    ("variant_filtration", "bcftools"): {"bcftools": {"exclude": "FILTER ~ 'low_depth'"}},
    ("variant_filtration", "regions"): {"regions": {"action": "exclude"}},
    ("variant_filtration", "vembrane"): {"vembrane": {"expressions": {"some_filter": "True"}}},
    ("variant_filtration", "varlociraptor"): {
        "varlociraptor": {"events": ["somatic_tumor_high", "somatic_tumor_medium"]},
    },
    ("wgs_cnv_export_external", None): VARIANT_EXPORT_EXTERNAL,
    ("wgs_sv_export_external", None): VARIANT_EXPORT_EXTERNAL,
}


def _external_file(
    name: str, data_type: str, tags: list[str], key: str, path: Fixture = REFERENCE
) -> dict[str, Any]:
    config = {"produces": {"type": data_type, "tags": tags}, "files": {key: path}}
    return {"step": "external_data", "name": name, "config": config}


def _external_files(name: str, data_type: str, tags: list[str], files: dict) -> dict[str, Any]:
    config = {"produces": {"type": data_type, "tags": tags}, "files": files}
    return {"step": "external_data", "name": name, "config": config}


#: The files of the mehari database fixture, as varfish_export read them below its directory
MEHARI_DB = {
    "hgnc_tsv": Fixture("mehari_db", "hgnc.tsv"),
    "transcript_db": Fixture("mehari_db", "grch37/seqvars/txs.bin.zst"),
    "clinvar_db": Fixture("mehari_db", "grch37/seqvars/clinvar/rocksdb"),
    "frequency_db": Fixture("mehari_db", "grch37/seqvars/frequencies/rocksdb"),
}


def _external_kits(name: str) -> dict[str, Any]:
    regions = {"targets": Fixture("cnvkit_targets"), "baits": Fixture("cnvkit_targets")}
    config = {"produces": {"type": "regions", "tags": ["kits"]}, "kits": {"MedExome": regions}}
    return {"step": "external_data", "name": name, "config": config}


def _external_vcfs(name: str, tags: list[str]) -> dict[str, Any]:
    pattern = {"vcf": r".+\.vcf\.gz", "vcf_tbi": r".+\.vcf\.gz\.tbi"}
    config = {"produces": {"type": "variants", "tags": tags}, "search_patterns": [pattern]}
    config["search_paths"] = Fixture("raw_data_dir")
    return {"step": "external_data", "name": name, "config": config}


#: external_data tasks besides the step's own one: the reference data, regions, population variants,
#: the regions of the kit of the test libraries (MedExome), and the inputs of the external export
#: steps. One task with several tags serves several kinds of regions; each task names the fixture
#: that the replaced config fields had, so that job inputs stay the same. features and regions
#: name the reference FASTA, as no test reads them; features also serves RNA-SeQC. dbsnp and
#: population_variants name a VCF, and mappable_regions an empty BED, whose indexes the jobs read
#: as well.
EXTRA_TASKS: list[dict[str, Any]] = [
    _external_file("genome", "raw", ["reference", "dna"], "fasta"),
    _external_file("features", "raw", ["features", "rnaseqc"], "gtf"),
    _external_file("dbsnp", "variants", ["dbsnp"], "vcf", VCF),
    _external_file("exons", "regions", ["exons"], "bed", PLACEHOLDER),
    _external_file("regions", "regions", ["microsatellites"], "bed"),
    _external_file("tagged_regions", "regions", ["par", "snp_loci", "genes"], "bed", PLACEHOLDER),
    _external_file("mappable_regions", "regions", ["mappable"], "bed", MAPPABLE),
    _external_file("population_variants", "variants", ["population", "common"], "vcf", VCF),
    _external_files(
        "salmon_index",
        "index",
        ["salmon", "rna"],
        {"index": Fixture("star_index"), "transcript_to_gene": PLACEHOLDER},
    ),
    _external_files("star_index", "index", ["star", "rna"], {"index": Fixture("star_index")}),
    _external_files("vep_cache", "database", ["vep"], {"cache": PLACEHOLDER}),
    _external_files("mehari_transcripts", "database", ["mehari"], {"transcript_db": [PLACEHOLDER]}),
    _external_files("mehari_db", "database", ["mehari"], MEHARI_DB),
    _external_files("arriba_database", "database", ["arriba"], {"blacklist": PLACEHOLDER}),
    _external_files(
        "star_fusion_database",
        "database",
        ["star_fusion"],
        {"ctat_resource_lib": Fixture("star_index")},
    ),
    _external_files(
        "control_freec_database", "database", ["control_freec"], {"chr_len": PLACEHOLDER}
    ),
    _external_files(
        "expansionhunter_database",
        "database",
        ["expansionhunter"],
        {"catalog": REFERENCE, "annotation": REFERENCE},
    ),
    _external_kits("kits"),
    _external_vcfs("external_vcf", ["germline", "snv", "indel"]),
    _external_vcfs("external_cnv", ["germline", "cnv"]),
    _external_vcfs("external_sv", ["germline", "sv"]),
    # The inputs of neoepitope prediction that other closures do not need
    _external_file("proteome", "raw", ["proteome"], "fasta", PLACEHOLDER),
    # FDR control of the germline varlociraptor calls per pedigree, and its calls report with IGV
    # views of the top calls, which read the BAMs and the genome
    {
        "step": "variant_filtration",
        "name": "variant_filtration_varlociraptor_germline",
        "config": {
            "depends_on": {"variants": "variant_calling_varlociraptor_germline"},
            "tool": "varlociraptor",
            "varlociraptor": {"events": ["present", "denovo_*", "recessive", "dominant"]},
        },
    },
    {
        "step": "variant_report",
        "name": "variant_report_igv",
        "config": {
            "depends_on": {
                "variants": "variant_filtration_varlociraptor_germline",
                "alignments": "ngs_mapping_bwa",
                "reference": "genome",
            },
            "tool": "datavzrd",
            "igv": {},
        },
    },
    # Candidate calling per pedigree, as the germline varlociraptor closure reads it
    {
        "step": "variant_calling",
        "name": "variant_calling_freebayes_germline",
        "config": {
            "depends_on": {"alignments": "ngs_mapping_bwa", "reference": "genome"},
            "group_by": "related",
            "tool": "freebayes",
            "freebayes": {"num_jobs": 2},
        },
    },
    {
        "step": "variant_calling",
        "name": "variant_calling_varlociraptor_germline",
        "config": {
            "depends_on": {
                "alignments": "ngs_mapping_bwa",
                "reference": "genome",
                "candidates": "variant_calling_freebayes_germline",
            },
            "tool": "varlociraptor",
            "varlociraptor": {"scenario": "germline", "num_jobs": 2},
        },
    },
    # Structural variants of the pedigree as candidates
    {
        "step": "variant_calling",
        "name": "variant_calling_varlociraptor_sv",
        "config": {
            "depends_on": {
                "alignments": "ngs_mapping_bwa",
                "reference": "genome",
                "candidates": "sv_calling_targeted_delly2",
            },
            "tool": "varlociraptor",
            "varlociraptor": {"scenario": "germline", "num_jobs": 2},
        },
    },
    {
        "step": "variant_annotation",
        "name": "variant_annotation_pvactools",
        "config": {
            "depends_on": {
                "variants": "variant_calling_mutect2",
                "reference": "genome",
                "database": "vep_cache",
            },
            "tool": "vep",
            "vep": {
                "pick": False,
                "plugins": [
                    {"name": "Frameshift", "path": PLACEHOLDER},
                    {"name": "Wildtype", "path": REFERENCE},
                ],
            },
        },
    },
    {
        "step": "variant_phasing",
        "name": "variant_phasing_pvactools",
        "config": {
            "depends_on": {
                "combined_variants": "combine_variants_annotated",
                "alignments": "ngs_mapping_bwa",
                "reference": "genome",
            },
            "tool": "gatk_read_backed_phasing_tumor",
        },
    },
    {
        "step": "combine_variants",
        "name": "combine_variants_annotated",
        "config": {
            "depends_on": {
                "somatic_variants": "variant_annotation_pvactools",
                "germline_variants": "variant_filtration_bcftools",
                "reference": "genome",
            }
        },
    },
    # The other sources of reference_download; the generated task downloads from Ensembl
    {
        "step": "reference_download",
        "name": "reference_download_refseq",
        "config": {
            "source": "refseq",
            "refseq": {"assembly_accession": "GCF_000001405.40", "assembly_name": "GRCh38.p14"},
        },
    },
    {
        "step": "reference_download",
        "name": "reference_download_ucsc",
        "config": {"source": "ucsc", "ucsc": {"db": "hg38"}},
    },
]


def _unset(value: Any) -> bool:
    return value is None or value in ("", "AUTO") or value in ([], ["AUTO"])


def _resolve(value: Any, fixtures: dict[str, Any]) -> Any:
    if isinstance(value, Fixture):
        path = copy.deepcopy(fixtures[value.name])
        return str(Path(path) / value.subpath) if value.subpath else path
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


def get_possible_variants(workflow_cls: type, step_name: str) -> list[str]:
    """Return the values of the field that splits a step into tasks; empty if it has none."""
    if step_name not in VARIANT_FIELDS:
        return []
    annotation = workflow_cls.config_model_class.model_fields[VARIANT_FIELDS[step_name]].annotation
    return [str(item.value) for item in annotation]


def tool_values(model: type) -> list[str]:
    """Return the values of a step model's ``tool`` field (enum or ``Literal``), if it has one."""
    tool_field = model.model_fields.get("tool")
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


#: The database task of each (step, tool); variant_annotation reads a list of transcript
#: databases, varfish_export the whole mehari database.
DATABASE_TASKS: dict[tuple[str, str | None], str] = {
    ("repeat_expansion", "expansionhunter"): "expansionhunter_database",
    ("somatic_gene_fusion_calling", "arriba"): "arriba_database",
    ("somatic_gene_fusion_calling", "star_fusion"): "star_fusion_database",
    ("somatic_wgs_cnv_calling", "control_freec"): "control_freec_database",
    ("varfish_export", None): "mehari_db",
    ("variant_annotation", "mehari"): "mehari_transcripts",
    ("variant_annotation", "vep"): "vep_cache",
}


def wire_dependencies(step_name: str, tool: str | None) -> dict[str, str | list[str]]:
    """Return the ``depends_on`` mapping of a generated task."""
    depends_on: dict[str, str | list[str]] = {}
    for field in dependency_fields(WORKFLOW_REGISTRY[step_name]):
        upstream = STEP_DEPENDENCY_TASKS.get((step_name, field), DEPENDENCY_TASKS.get(field))
        if upstream:
            depends_on[field] = upstream

    # featurecounts and dupradar read the strandedness decision of the strandedness task.
    if step_name == "gene_expression_quantification" and tool in ("featurecounts", "dupradar"):
        depends_on["strandedness"] = "gene_expression_quantification_strandedness"

    # arcasHLA types RNA, which only the STAR mapping provides.
    if step_name == "hla_typing" and tool == "arcashla":
        depends_on["alignments"] = "ngs_mapping_star"

    # A mapper reads the index of its own tool.
    if step_name == "ngs_mapping" and tool is not None:
        depends_on["index"] = f"reference_index_{tool}"

    # The purecn panel of normals builds on the mutect2 one, the cnvkit one uses the access regions.
    # The cnvkit panel stays whole genome, as somatic_wgs_cnv_calling needs it.
    if step_name == "panel_of_normals" and tool == "cnvkit":
        depends_on.pop("kits", None)
    if step_name == "panel_of_normals" and tool == "purecn":
        depends_on["panel_of_normals"] = f"{step_name}_mutect2"
    if step_name == "panel_of_normals" and tool == "cnvkit":
        depends_on["access"] = f"{step_name}_access"

    # Only the regions filter reads regions.
    if step_name == "variant_filtration" and tool == "regions":
        depends_on["regions"] = "regions"

    # FDR control reads varlociraptor calls.
    if step_name == "variant_filtration" and tool == "varlociraptor":
        depends_on["variants"] = "variant_calling_varlociraptor"

    # A tool reads its own database.
    database = DATABASE_TASKS.get((step_name, tool))
    if database and "database" in dependency_fields(WORKFLOW_REGISTRY[step_name]):
        depends_on["database"] = database

    # Read-backed phasing refines the result of phasing by transmission.
    if step_name == "variant_phasing" and tool == "gatk_read_backed_phasing":
        depends_on["variants"] = "variant_phasing_gatk_phase_by_transmission"

    # Mutect2 calling filters with the mutect2 panel of normals.
    if step_name == "variant_calling" and tool == "mutect2":
        depends_on["panel_of_normals"] = "panel_of_normals_mutect2"

    # The gCNV callers read the model of their own kind of intervals.
    if step_name == "sv_calling_targeted" and tool == "gcnv":
        depends_on["gcnv_model"] = "gcnv_model_targeted"
    if step_name == "sv_calling_wgs" and tool == "gcnv":
        depends_on["gcnv_model"] = "gcnv_model_wgs"

    # cnvkit and purecn read the panel of normals built with the same tool.
    if step_name == "somatic_targeted_seq_cnv_calling" and tool in ("cnvkit", "purecn"):
        depends_on["panel_of_normals"] = f"panel_of_normals_{tool}"
    if step_name == "somatic_wgs_cnv_calling" and tool == "cnvkit":
        depends_on["panel_of_normals"] = "panel_of_normals_cnvkit"
    return depends_on


def build_all_tasks(base_config: dict[str, Any], base_config_path: Path) -> list[dict[str, Any]]:
    """Return one task per step and tool, with filled config and wired dependencies."""
    fixtures = fixture_paths(base_config, base_config_path)
    tasks: list[dict[str, Any]] = []
    for step_name, cls in sorted(WORKFLOW_REGISTRY.items()):
        variants = get_possible_variants(cls, step_name)
        for tool in variants or tool_values(cls.config_model_class) or [None]:
            config = parse_default_step_config(cls)
            if variants:
                config[VARIANT_FIELDS[step_name]] = tool
            elif tool is not None:
                config["tool"] = tool
                # The selected tool's section must be present, even if empty.
                if tool in cls.config_model_class.model_fields:
                    config.setdefault(tool, {})
            fill_config(config, TASK_CONFIG.get((step_name, None), {}), fixtures)
            if tool is not None:
                fill_config(config, TASK_CONFIG.get((step_name, tool), {}), fixtures)
            if depends_on := wire_dependencies(step_name, tool):
                config["depends_on"] = depends_on

            name = f"{step_name}_{tool}" if tool else step_name
            config = validate_step_config(name, cls, config)
            if depends_on:
                # Validation dumps every depends_on field; keep the wired ones that are read.
                read = config.get("depends_on", {})
                config["depends_on"] = {k: v for k, v in depends_on.items() if read.get(k)}
            tasks.append({"step": step_name, "name": name, "config": config})

    tasks += [_resolve(extra, fixtures) for extra in EXTRA_TASKS]
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
    config = {
        "project": copy.deepcopy(base_config["project"]),
        "data_sets": copy.deepcopy(base_config.get("data_sets", {})),
        "tasks": tasks,
    }
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
