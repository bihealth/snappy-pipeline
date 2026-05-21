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
from typing import Any, get_args, get_origin

import pydantic
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import DataSignature

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
            validated = config_model(**candidate)
            return validated.model_dump(exclude_none=True), notes
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
    dep_instance = dep_model()
    dep_dict = dep_instance.model_dump()
    return {k: v for k, v in dep_dict.items()}


def _guess_bwa_index_from_reference(base_config: dict[str, Any]) -> str:
    static_data = base_config.get("static_data_config", {})
    if isinstance(static_data, dict):
        ref = static_data.get("reference", {})
        if isinstance(ref, dict):
            path = ref.get("path")
            if isinstance(path, str) and path:
                return path
    return "AUTO"


def bootstrap_step_config(
    step_name: str,
    step_config: dict[str, Any],
    base_config: dict[str, Any],
) -> dict[str, Any]:
    cfg = copy.deepcopy(step_config)

    if step_name == "ngs_mapping":
        cfg["tool"] = "bwa"
        cfg.setdefault(
            "target_coverage_report",
            {"enabled": False, "path_target_interval_list_mapping": []},
        )
        cfg.setdefault("bwa", {})
        if isinstance(cfg["bwa"], dict):
            cfg["bwa"].setdefault("path_index", _guess_bwa_index_from_reference(base_config))

    if step_name == "variant_calling":
        cfg["tool"] = "bcftools_call"
        cfg.setdefault("bcftools_call", {})
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
) -> str | None:
    candidates: list[str] = []
    for step_name, cls in workflow_items:
        if step_name == consumer_step:
            continue
        produces = getattr(cls, "produces", []) or []
        if any(sig.satisfies(requirement) for sig in produces):
            candidates.append(step_name)

    if not candidates:
        return None

    candidates = sorted(candidates, key=lambda c: _candidate_score(c, requirement), reverse=True)
    return candidates[0]


def build_all_tasks(base_config: dict[str, Any], base_config_path: Path) -> list[dict[str, Any]]:
    workflow_items = sorted(WORKFLOW_REGISTRY.items())
    all_steps = {name for name, _ in workflow_items}

    generation_notes: list[str] = []
    tasks: list[dict[str, Any]] = []
    for step_name, cls in workflow_items:
        step_config = parse_default_step_config(cls, step_name)
        step_config = bootstrap_step_config(step_name, step_config, base_config)
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
                "name": step_name,
                "config": step_config,
            }
        )

    tasks_by_step = {t["step"]: t for t in tasks}

    for step_name, cls in workflow_items:
        task = tasks_by_step[step_name]
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
                depends_on[logical_name] = default_target
            elif logical_name in all_steps and logical_name != step_name:
                depends_on[logical_name] = logical_name
            elif logical_name in LOGICAL_DEP_ALIASES:
                alias_target = LOGICAL_DEP_ALIASES[logical_name]
                if alias_target in all_steps and alias_target != step_name:
                    depends_on[logical_name] = alias_target

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
            producer = find_producer_for_requirement(selected_req, workflow_items, step_name)
            if producer and producer != step_name:
                depends_on[logical_name] = producer

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
    normalize_data_set_paths(config, base_config_path)
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=Path(
            ".tests/test-workflow/pipelines/snappy-cancer_wes/.snappy_pipeline/config.yaml"
        ),
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
