#!/usr/bin/env python3
"""Generate task-based config files that include all registered workflows.

The generator builds one task per workflow step, derives each task's config from
`default_config_yaml()`, and wires `depends_on` mappings systematically using:
1) typed `depends_on` defaults from each step config model,
2) explicit logical-name aliases for non-step dependency names, and
3) consumes/produces signature matching as a fallback.
"""

from __future__ import annotations

import argparse
import copy
import enum
from pathlib import Path
import typing
from typing import Any, get_args, get_origin

import pydantic
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import DataSignature, ExpectedPathSchema

yaml = ruamel_yaml.YAML()
yaml.default_flow_style = False
yaml.indent(sequence=4, offset=2)


# Logical dependency names that are not step names but commonly appear in depends_on models.
LOGICAL_DEP_ALIASES: dict[str, str] = {
    "somatic_variant": "somatic_variant_calling",
    "somatic_variants": "somatic_variant_calling",
    "cnv_calling": "somatic_wgs_cnv_calling",
    "copy_number": "somatic_wgs_cnv_calling",
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


def infer_link_in_path(base_config: dict[str, Any], base_config_path: Path) -> str:
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


def get_dep_defaults(workflow_cls: type) -> dict[str, str | None]:
    config_model = getattr(workflow_cls, "config_model_class", None)
    if not config_model:
        return {}

    model_fields = getattr(config_model, "model_fields", {})
    dep_field = model_fields.get("depends_on")
    if dep_field is None:
        return {}

    dep_model = dep_field.annotation
    dep_dict = {}
    from pydantic_core import PydanticUndefined

    for name, field in getattr(dep_model, "model_fields", {}).items():
        val = field.default
        if val is PydanticUndefined:
            dep_dict[name] = ""
        else:
            dep_dict[name] = val
    return dep_dict


def _guess_bwa_index_from_reference(base_config: dict[str, Any]) -> str:
    static_data = base_config.get("static_data_config", {})
    if isinstance(static_data, dict):
        ref = static_data.get("reference", {})
        if isinstance(ref, dict):
            path = ref.get("path")
            if isinstance(path, str) and path:
                return path
    return "AUTO"


def _existing_placeholder_file() -> str:
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "test.fai"
    if candidate.exists():
        return str(candidate)
    return str(__file__)


def _star_index_fixture_dir() -> str:
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "snappy_pipeline" / "fixtures" / "star_index"
    if candidate.exists():
        return str(candidate)
    return str(repo_root)


def _cnvkit_targets_bed() -> str:
    """Return path to cnvkit target regions fixture BED file."""
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "snappy_pipeline" / "fixtures" / "cnvkit_targets.bed"
    if candidate.exists():
        return str(candidate)
    # Also check test fixture directory
    candidate2 = repo_root / "tests" / "snappy_pipeline" / "fixtures" / "cnvkit_targets.bed"
    if candidate2.exists():
        return str(candidate2)
    return str(__file__)


def _cnvkit_antitargets_bed() -> str:
    """Return path to cnvkit antitarget regions fixture BED file."""
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "snappy_pipeline" / "fixtures" / "cnvkit_antitargets.bed"
    if candidate.exists():
        return str(candidate)
    # Also check test fixture directory
    candidate2 = repo_root / "tests" / "snappy_pipeline" / "fixtures" / "cnvkit_antitargets.bed"
    if candidate2.exists():
        return str(candidate2)
    return str(__file__)


def _gcnv_ploidy_model_dir() -> str:
    """Return path to gCNV ploidy model fixture."""
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "snappy_pipeline" / "fixtures" / "gcnv_models" / "ploidy_model"
    if candidate.exists():
        return str(candidate)
    return str(repo_root)


def _gcnv_call_model_pattern() -> str:
    """Return path pattern for gCNV call model fixtures."""
    repo_root = Path(__file__).resolve().parent.parent
    candidate = repo_root / "snappy_pipeline" / "fixtures" / "gcnv_models" / "call_model_*"
    # Return the base directory, the pattern will be expanded by get_model_dir_list
    return str(candidate)


def _get_gcnv_precomputed_models(workflow_type: str = "targeted") -> list[dict[str, str]]:
    """Generate precomputed model paths for gCNV testing.

    Args:
        workflow_type: Either "targeted" (uses "default" library name) or "wgs" (uses "wgs" library name)
    """
    library_name = "wgs" if workflow_type == "wgs" else "default"
    return [
        {
            "library": library_name,
            "contig_ploidy": _gcnv_ploidy_model_dir(),
            "model_pattern": _gcnv_call_model_pattern(),
        }
    ]


def _guess_reference_from_static_data(base_config: dict[str, Any]) -> str:
    static_data = base_config.get("static_data_config", {})
    if isinstance(static_data, dict):
        ref = static_data.get("reference", {})
        if isinstance(ref, dict):
            path = ref.get("path")
            if isinstance(path, str) and path:
                return path
    return _existing_placeholder_file()


def _resolve_path(path: str, base_config_path: Path | None) -> str:
    """Resolve a path to absolute if it's relative and base_config_path is provided."""
    if not path or path.startswith("/") or path.startswith("AUTO"):
        return path
    if base_config_path is None:
        return path
    return str((base_config_path.parent / path).resolve())


def bootstrap_step_config(
    step_name: str,
    step_config: dict[str, Any],
    base_config: dict[str, Any],
    base_config_path: Path | None = None,
) -> dict[str, Any]:
    cfg = copy.deepcopy(step_config)

    if step_name == "reference_index":
        tool = cfg.get("tool") or "bwa"
        cfg["tool"] = tool
        if tool == "star":
            cfg["reference_molecule"] = "rna"

    if step_name == "ngs_mapping":
        tool = cfg.get("tool") or "bwa"
        cfg["tool"] = tool
        cfg.setdefault(
            "target_coverage_report",
            {"enabled": False, "path_target_interval_list_mapping": []},
        )
        cfg.setdefault(tool, {})
        if tool in ("bwa", "bwa_mem2", "minimap2") and isinstance(cfg[tool], dict):
            cfg[tool].setdefault("path_index", _guess_bwa_index_from_reference(base_config))
        elif tool == "star" and isinstance(cfg["star"], dict):
            cfg["star"].setdefault("path_index", _star_index_fixture_dir())
            cfg.setdefault("strandedness", {})
            if isinstance(cfg["strandedness"], dict):
                ref_path = _guess_reference_from_static_data(base_config).replace(
                    ".fa", ".exon.bed"
                )
                cfg["strandedness"].setdefault(
                    "path_exon_bed",
                    _resolve_path(ref_path, base_config_path),
                )
                cfg["strandedness"].setdefault("strand", -1)
                cfg["strandedness"].setdefault("threshold", 0.85)
        elif tool == "mbcs":
            if isinstance(cfg["mbcs"], dict):
                cfg["mbcs"].setdefault("mapping_tool", "bwa")
            cfg.setdefault("bwa", {})
            if isinstance(cfg["bwa"], dict):
                cfg["bwa"].setdefault("path_index", _guess_bwa_index_from_reference(base_config))
            cfg.setdefault("bqsr", {})
            if isinstance(cfg["bqsr"], dict):
                ref_path = _guess_reference_from_static_data(base_config)
                cfg["bqsr"].setdefault("common_variants", _resolve_path(ref_path, base_config_path))

    if step_name == "ngs_data_qc":
        tool = cfg.get("tool") or "fastqc"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "picard" and isinstance(cfg["picard"], dict):
            programs = cfg["picard"].get("programs")
            if not isinstance(programs, list) or not programs:
                cfg["picard"]["programs"] = ["CollectAlignmentSummaryMetrics"]

    if step_name == "somatic_targeted_seq_cnv_calling":
        # HRD requires sequenza; sequenza also produces _dnacopy.seg used by cnv_checking
        tool = cfg.get("tool") or "sequenza"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "cnvkit" and isinstance(cfg.get("cnvkit"), dict):
            # Use fixture BED files for target/antitarget regions
            if cfg["cnvkit"].get("path_target") in (None, "", "AUTO"):
                cfg["cnvkit"]["path_target"] = _cnvkit_targets_bed()
            if cfg["cnvkit"].get("path_antitarget") in (None, "", "AUTO"):
                cfg["cnvkit"]["path_antitarget"] = _cnvkit_antitargets_bed()
            # path_panel_of_normals is no longer a config field; it comes from depends_on.panel_of_normals
        elif tool == "purecn" and isinstance(cfg.get("purecn"), dict):
            # path_panel_of_normals / path_intervals / path_mapping_bias are no longer config fields;
            # they come from depends_on.panel_of_normals resolved at runtime.
            if not isinstance(cfg["purecn"].get("path_container"), str) or cfg["purecn"].get(
                "path_container"
            ) in ("", "AUTO"):
                cfg["purecn"]["path_container"] = _existing_placeholder_file()

    if step_name == "somatic_wgs_cnv_calling":
        # cnvkit produces _dnacopy.seg expected by somatic_cnv_checking
        tool = cfg.get("tool") or "cnvkit"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})

    if step_name == "sv_calling_targeted":
        tool = cfg.get("tool") or "delly2"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "gcnv" and isinstance(cfg.get("gcnv"), dict):
            cfg["gcnv"].setdefault(
                "precomputed_model_paths", _get_gcnv_precomputed_models("targeted")
            )

    if step_name == "sv_calling_wgs":
        tool = cfg.get("tool") or "delly2"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "gcnv" and isinstance(cfg.get("gcnv"), dict):
            cfg["gcnv"].setdefault("precomputed_model_paths", _get_gcnv_precomputed_models("wgs"))

    if step_name == "repeat_expansion":
        placeholder = _guess_reference_from_static_data(base_config)
        placeholder = _resolve_path(placeholder, base_config_path)
        if not isinstance(cfg.get("repeat_catalog"), str) or cfg.get("repeat_catalog") in (
            "",
            "AUTO",
        ):
            cfg["repeat_catalog"] = placeholder
        if not isinstance(cfg.get("repeat_annotation"), str) or cfg.get("repeat_annotation") in (
            "",
            "AUTO",
        ):
            cfg["repeat_annotation"] = placeholder

    if step_name == "panel_of_normals":
        tool = cfg.get("tool") or "mutect2"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "mutect2" and isinstance(cfg["mutect2"], dict):
            if not isinstance(cfg["mutect2"].get("germline_resource"), str) or cfg["mutect2"].get(
                "germline_resource"
            ) in ("", "AUTO"):
                ref_path = _guess_reference_from_static_data(base_config)
                cfg["mutect2"]["germline_resource"] = _resolve_path(ref_path, base_config_path)
        elif tool == "cnvkit" and isinstance(cfg["cnvkit"], dict):
            # Empty target path puts CNVkit into WGS mode and avoids unresolved external target BEDs.
            cfg["cnvkit"]["path_target"] = ""
        elif tool == "purecn" and isinstance(cfg["purecn"], dict):
            if not isinstance(cfg["purecn"].get("path_bait_regions"), str) or cfg["purecn"].get(
                "path_bait_regions"
            ) in ("", "AUTO"):
                ref_path = _guess_reference_from_static_data(base_config)
                cfg["purecn"]["path_bait_regions"] = _resolve_path(ref_path, base_config_path)
            cfg["purecn"]["path_normals_list"] = ""
            # path_genomicsDB removed: the genomicsDB is now a tracked Snakemake input derived
            # from depends_on.panel_of_normals, not a bare config path.

    if step_name == "somatic_gene_fusion_calling":
        tool = cfg.get("tool") or "arriba"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "arriba" and isinstance(cfg["arriba"], dict):
            cfg["arriba"].setdefault("path_index", _star_index_fixture_dir())

    if step_name in ("somatic_variant_annotation", "variant_annotation"):
        tool = cfg.get("tool") or "vep"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "mehari" and isinstance(cfg["mehari"], dict):
            if not isinstance(cfg["mehari"].get("reference"), str) or cfg["mehari"].get(
                "reference"
            ) in ("", "AUTO"):
                ref_path = _guess_reference_from_static_data(base_config)
                cfg["mehari"]["reference"] = _resolve_path(ref_path, base_config_path)
            if not isinstance(cfg["mehari"].get("transcripts"), list):
                cfg["mehari"]["transcripts"] = [_existing_placeholder_file()]
            elif cfg["mehari"].get("transcripts") == ["AUTO"]:
                cfg["mehari"]["transcripts"] = [_existing_placeholder_file()]

    if step_name == "somatic_msi_calling":
        tool = cfg.get("tool") or "mantis_msi2"
        cfg["tool"] = tool
        if not isinstance(cfg.get("loci_bed"), str) or cfg.get("loci_bed") in ("", "AUTO"):
            ref_path = _guess_reference_from_static_data(base_config)
            cfg["loci_bed"] = _resolve_path(ref_path, base_config_path)

    if step_name == "targeted_seq_mei_calling":
        tool = cfg.get("tool") or "scramble"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        if tool == "scramble" and isinstance(cfg["scramble"], dict):
            if not isinstance(cfg["scramble"].get("blast_ref"), str) or cfg["scramble"].get(
                "blast_ref"
            ) in ("", "AUTO"):
                ref_path = _guess_reference_from_static_data(base_config)
                cfg["scramble"]["blast_ref"] = _resolve_path(ref_path, base_config_path)

    if step_name in (
        "variant_export_external",
        "wgs_cnv_export_external",
        "wgs_sv_export_external",
    ):
        # These require list-typed non-empty search config and existing DB/serializer files.
        if not isinstance(cfg.get("search_paths"), list) or not cfg.get("search_paths"):
            cfg["search_paths"] = ["/tmp"]
        if not isinstance(cfg.get("search_patterns"), list) or not cfg.get("search_patterns"):
            cfg["search_patterns"] = [{"vcf": "*.vcf.gz"}]
        placeholder_file = _existing_placeholder_file()
        for key in ("path_refseq_ser", "path_ensembl_ser", "path_db"):
            value = cfg.get(key)
            if not isinstance(value, str) or not value or value == "AUTO":
                cfg[key] = placeholder_file

    if step_name == "variant_calling":
        tool = cfg.get("tool") or "bcftools_call"
        cfg["tool"] = tool
        cfg.setdefault(tool, {})
        cfg.setdefault("baf_file_generation", {"enabled": False, "min_dp": 10})
        cfg.setdefault("bcftools_stats", {"enabled": False})
        cfg.setdefault("jannovar_stats", {"enabled": False, "path_ser": "AUTO"})
        cfg.setdefault(
            "bcftools_roh",
            {
                "enabled": False,
                "path_af_file": "AUTO",
                "path_targets": None,
                "ignore_homref": False,
                "skip_indels": False,
                "rec_rate": 1e-8,
            },
        )

    if step_name == "variant_filtration":
        tool = cfg.get("tool") or "bcftools"
        cfg["tool"] = tool
        if tool == "bcftools":
            cfg.setdefault("bcftools", {})
            if isinstance(cfg["bcftools"], dict):
                cfg["bcftools"].setdefault("exclude", "FILTER ~ 'low_depth'")
        elif tool == "regions":
            cfg.setdefault("regions", {})
            if isinstance(cfg["regions"], dict):
                cfg["regions"].setdefault("exclude", "FILTER ~ 'low_depth'")
        elif tool == "vembrane":
            cfg.setdefault("vembrane", {})
            if isinstance(cfg["vembrane"], dict):
                cfg["vembrane"].setdefault("expressions", {"some_filter": "True"})

    return cfg


def ensure_explicit_selected_tool_config(
    step_name: str,
    workflow_cls: type,
    step_config: dict[str, Any],
) -> tuple[dict[str, Any], list[str]]:
    config_model = getattr(workflow_cls, "config_model_class", None)
    if not config_model:
        return step_config, []

    model_fields = getattr(config_model, "model_fields", {})
    if "tool" not in model_fields:
        return step_config, []

    tool_field = model_fields.get("tool")
    selected_tool = step_config.get("tool")
    if selected_tool is None and tool_field is not None:
        selected_tool = tool_field.default
    if isinstance(selected_tool, enum.Enum):
        selected_tool = selected_tool.value
    if not isinstance(selected_tool, str) or not selected_tool:
        return step_config, []

    # Most workflow models name the tool-specific section after the tool value.
    if selected_tool not in model_fields:
        return step_config, []

    cfg = copy.deepcopy(step_config)
    if "tool" not in cfg:
        cfg["tool"] = selected_tool
        note = f"{step_name}: added explicit tool selection tool: {selected_tool}"
    else:
        note = None

    if selected_tool in cfg:
        return cfg, ([note] if note else [])

    cfg[selected_tool] = {}
    notes = [f"{step_name}: added explicit selected-tool section {selected_tool}: {{}}"]
    if note:
        notes.insert(0, note)
    return cfg, notes


def _candidate_score(candidate_step: str, requirement: DataSignature) -> tuple[int, int, str]:
    tags = getattr(requirement, "tags", frozenset())
    positive_tags = [t for t in tags if isinstance(t, str) and not t.startswith("-")]
    score = len(positive_tags)

    # Mild domain preference: for DNA alignments choose mapping first.
    if requirement.type.value == "alignments" and candidate_step == "ngs_mapping":
        score += 2

    return (score, len(candidate_step), candidate_step)


def find_producer_for_requirement(
    requirement: DataSignature,
    workflow_items: list[tuple[str, type]],
    consumer_step: str,
    expected_schema: type[pydantic.BaseModel] | None = None,
) -> str | None:
    candidates: list[str] = []
    for step_name, cls in workflow_items:
        if step_name == consumer_step:
            continue
        produces = getattr(cls, "produces", []) or []
        if any(sig.satisfies(requirement) for sig in produces):
            if expected_schema is not None:
                try:
                    out_paths = cls.get_output_paths(signature=requirement)
                    if isinstance(out_paths, dict):
                        if not all(
                            field_name in out_paths for field_name in expected_schema.model_fields
                        ):
                            continue
                except Exception:
                    continue
            candidates.append(step_name)

    if not candidates:
        return None

    candidates = sorted(candidates, key=lambda c: _candidate_score(c, requirement), reverse=True)
    return candidates[0]


def get_possible_tools(workflow_cls: type) -> list[str]:
    config_model = getattr(workflow_cls, "config_model_class", None)
    if not config_model:
        return []
    model_fields = getattr(config_model, "model_fields", {})
    tool_field = model_fields.get("tool")
    if tool_field is None:
        return []

    ann = tool_field.annotation
    # Unwrap Annotated, Union, Optional
    while True:
        origin = get_origin(ann)
        args = get_args(ann)
        if origin is not None and getattr(origin, "__name__", None) == "Annotated":
            ann = args[0]
            continue
        if origin is typing.Union or (
            origin is not None and getattr(origin, "__name__", None) in ("Union", "UnionType")
        ):
            non_none_args = [a for a in args if a is not type(None)]
            if non_none_args:
                ann = non_none_args[0]
                continue
        break

    # Now check if ann is Enum or Literal
    origin = get_origin(ann)
    if origin is typing.Literal or getattr(origin, "__name__", None) == "Literal":
        return [str(val) for val in get_args(ann)]

    if isinstance(ann, type) and issubclass(ann, enum.Enum):
        return [str(item.value) for item in ann]

    return []


PREFERRED_DEFAULT_TOOLS: dict[str, str] = {
    # Prefer configs that can dryrun without heavy precomputed models.
    "somatic_targeted_seq_cnv_calling": "sequenza",
    "somatic_wgs_cnv_calling": "cnvkit",
    "sv_calling_targeted": "delly2",
}


def get_default_tool(step_name: str, workflow_cls: type) -> str | None:
    preferred = PREFERRED_DEFAULT_TOOLS.get(step_name)
    if preferred:
        possible = get_possible_tools(workflow_cls)
        if preferred in possible:
            return preferred

    config_model = getattr(workflow_cls, "config_model_class", None)
    if not config_model:
        return None
    model_fields = getattr(config_model, "model_fields", {})
    tool_field = model_fields.get("tool")
    if tool_field is None:
        return None
    default_val = tool_field.default
    from pydantic_core import PydanticUndefined

    if default_val is PydanticUndefined:
        default_val = None
    if default_val is not None:
        if isinstance(default_val, enum.Enum):
            return str(default_val.value)
        if isinstance(default_val, str):
            return default_val
    possible = get_possible_tools(workflow_cls)
    if possible:
        return possible[0]
    return None


def get_default_task_name(step_name: str, workflow_cls: type) -> str:
    """Return the task name for the default tool of a step.

    Uses the actual tool name (e.g., 'mutect2', 'bwa') instead of the generic '_default' suffix.
    This makes task names explicit and easier to debug.
    """
    default_tool = get_default_tool(step_name, workflow_cls)
    if default_tool:
        return f"{step_name}_{default_tool}"
    return step_name


def build_all_tasks(base_config: dict[str, Any], base_config_path: Path) -> list[dict[str, Any]]:
    workflow_items = sorted(WORKFLOW_REGISTRY.items())
    all_steps = {name for name, _ in workflow_items}

    step_to_default_task: dict[str, str] = {}
    for name, cls in workflow_items:
        step_to_default_task[name] = get_default_task_name(name, cls)

    generation_notes: list[str] = []
    tasks: list[dict[str, Any]] = []
    for step_name, cls in workflow_items:
        tools = get_possible_tools(cls)

        if tools:
            for tool_name in tools:
                task_name = f"{step_name}_{tool_name}"

                step_config = parse_default_step_config(cls, step_name)
                step_config["tool"] = tool_name
                step_config = bootstrap_step_config(
                    step_name, step_config, base_config, base_config_path
                )
                if step_name == "link_in" and "path" not in step_config:
                    step_config["path"] = infer_link_in_path(base_config, base_config_path)

                step_config, notes = validate_and_autofill_step_config(step_name, cls, step_config)
                generation_notes.extend(notes)
                step_config, notes = ensure_explicit_selected_tool_config(
                    step_name, cls, step_config
                )
                generation_notes.extend(notes)
                if notes:
                    step_config, notes = validate_and_autofill_step_config(
                        step_name, cls, step_config
                    )
                    generation_notes.extend(notes)

                tasks.append(
                    {
                        "step": step_name,
                        "name": task_name,
                        "config": step_config,
                    }
                )
        else:
            task_name = step_name
            step_config = parse_default_step_config(cls, step_name)
            step_config = bootstrap_step_config(
                step_name, step_config, base_config, base_config_path
            )
            if step_name == "link_in" and "path" not in step_config:
                step_config["path"] = infer_link_in_path(base_config, base_config_path)

            step_config, notes = validate_and_autofill_step_config(step_name, cls, step_config)
            generation_notes.extend(notes)
            step_config, notes = ensure_explicit_selected_tool_config(step_name, cls, step_config)
            generation_notes.extend(notes)
            if notes:
                step_config, notes = validate_and_autofill_step_config(step_name, cls, step_config)
                generation_notes.extend(notes)

            tasks.append(
                {
                    "step": step_name,
                    "name": task_name,
                    "config": step_config,
                }
            )

    for task in tasks:
        step_name = task["step"]
        cls = WORKFLOW_REGISTRY[step_name]
        depends_on: dict[str, str] = {}

        # 1) Use typed depends_on defaults first (if available).
        dep_defaults = get_dep_defaults(cls)
        for logical_name, default_target in dep_defaults.items():
            if (
                isinstance(default_target, str)
                and default_target
                and default_target in all_steps
                and default_target != step_name
            ):
                if (
                    step_name == "homologous_recombination_deficiency"
                    and default_target == "somatic_targeted_seq_cnv_calling"
                ):
                    depends_on[logical_name] = "somatic_targeted_seq_cnv_calling_sequenza"
                else:
                    depends_on[logical_name] = step_to_default_task[default_target]
            elif logical_name in all_steps and logical_name != step_name:
                if (
                    step_name == "homologous_recombination_deficiency"
                    and logical_name == "somatic_targeted_seq_cnv_calling"
                ):
                    depends_on[logical_name] = "somatic_targeted_seq_cnv_calling_sequenza"
                else:
                    depends_on[logical_name] = step_to_default_task[logical_name]
            elif logical_name in LOGICAL_DEP_ALIASES:
                alias_target = LOGICAL_DEP_ALIASES[logical_name]
                if alias_target in all_steps and alias_target != step_name:
                    if (
                        step_name == "homologous_recombination_deficiency"
                        and alias_target == "somatic_targeted_seq_cnv_calling"
                    ):
                        depends_on[logical_name] = "somatic_targeted_seq_cnv_calling_sequenza"
                    else:
                        depends_on[logical_name] = step_to_default_task[alias_target]

        # 2) Fill remaining typed dependency keys by consumes/produces matching.
        unresolved = [k for k in dep_defaults if k not in depends_on]
        required_reqs = [
            req for req, required in (getattr(cls, "consumes", {}) or {}).items() if required
        ]
        for logical_name in unresolved:
            # Prefer a requirement whose type name resembles the logical dependency name.
            selected_req = next(
                (
                    req
                    for req in required_reqs
                    if req.type.value in logical_name or logical_name in req.type.value
                ),
                None,
            )
            if selected_req is None and required_reqs:
                selected_req = required_reqs[0]
            if selected_req is None:
                continue
            expected_schema = None
            config_model = getattr(cls, "config_model_class", None)
            if config_model:
                dep_field = config_model.model_fields.get("depends_on")
                if dep_field is not None:
                    dep_model = dep_field.annotation
                    field_info = dep_model.model_fields.get(logical_name)
                    if field_info is not None:
                        for item in getattr(field_info, "metadata", []):
                            if isinstance(item, ExpectedPathSchema):
                                expected_schema = item.schema
                            elif isinstance(item, type) and issubclass(item, pydantic.BaseModel):
                                expected_schema = item

            producer = find_producer_for_requirement(
                selected_req, workflow_items, step_name, expected_schema=expected_schema
            )
            if producer and producer != step_name:
                depends_on[logical_name] = step_to_default_task[producer]

        # 3) Special handling for panel_of_normals:
        # - For purecn: depends_on.panel_of_normals must point to the mutect2 variant
        # - For other tools: remove any panel_of_normals dependency (it's only for purecn)
        task_config = task.get("config", {})
        if step_name == "panel_of_normals" and isinstance(task_config, dict):
            if task_config.get("tool") == "purecn":
                depends_on["panel_of_normals"] = f"{step_name}_mutect2"
            else:
                # Remove panel_of_normals dependency for non-purecn tools
                depends_on.pop("panel_of_normals", None)

        # 4) Special handling for somatic_targeted_seq_cnv_calling:
        # - For cnvkit: depends_on.panel_of_normals -> panel_of_normals_cnvkit
        # - For purecn: depends_on.panel_of_normals -> panel_of_normals_purecn
        # - For sequenza: no panel_of_normals dependency needed
        if step_name == "somatic_targeted_seq_cnv_calling" and isinstance(task_config, dict):
            tool_val = task_config.get("tool")
            if tool_val in ("cnvkit", "purecn"):
                depends_on["panel_of_normals"] = f"panel_of_normals_{tool_val}"
            else:
                depends_on.pop("panel_of_normals", None)

        if depends_on:
            task_config = task.get("config", {})
            if isinstance(task_config, dict):
                task_config["depends_on"] = depends_on

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
    config = {
        "static_data_config": copy.deepcopy(base_config.get("static_data_config", {})),
        "tasks": tasks,
        "data_sets": copy.deepcopy(base_config.get("data_sets", {})),
    }

    # Enrich static_data_config with placeholders for optional but frequently accessed fields
    static_data = config.get("static_data_config", {})
    if isinstance(static_data, dict):
        # Resolve any relative paths to absolute relative to the base config file
        for k, v in static_data.items():
            if isinstance(v, dict) and "path" in v and isinstance(v["path"], str):
                p = v["path"]
                if p and not p.startswith("/") and not p.startswith("AUTO"):
                    v["path"] = str((base_config_path.parent / p).resolve())

        ref_path = "AUTO"
        ref_obj = static_data.get("reference")
        if isinstance(ref_obj, dict) and ref_obj.get("path"):
            ref_path = ref_obj["path"]
        else:
            ref_path = _existing_placeholder_file()

        for key in ("cosmic", "dbsnp", "dbnsfp", "features"):
            if key not in static_data or static_data[key] is None:
                static_data[key] = {"path": ref_path}

    normalize_data_set_paths(config, base_config_path)
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=Path("tests/snappy_pipeline/fixtures/base_config.yaml"),
        help="Path to an existing config.yaml used as source for static_data_config and data_sets.",
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
