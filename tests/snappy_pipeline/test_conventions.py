# -*- coding: utf-8 -*-
"""Checks for the coding conventions in ``docs/dev_conventions.rst``."""

from __future__ import annotations

import ast
import enum
import importlib
import json
import os
import re
from collections.abc import Iterator
from pathlib import Path

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = REPO / "snappy_pipeline" / "workflows"
SNAPSHOTS = Path(__file__).resolve().parent / "snapshots"

#: Path convention violations that already exist (plans.md C3). The test fails for new ones and
#: for fixed ones, so the list only shrinks.
PATH_BASELINE = SNAPSHOTS / "path_convention_violations.json"

#: Input functions that a Snakefile calls with Snakemake's ``checkpoints`` object.
CHECKPOINT_INPUT_FUNCTIONS = {
    ("BuildGcnvTargetSeqModelStepPart", "_get_input_files_post_germline_calls"),
    ("BuildGcnvWgsModelStepPart", "_get_input_files_post_germline_calls"),
}


def _class_members() -> Iterator[tuple[Path, ast.ClassDef, ast.stmt]]:
    for path in sorted(WORKFLOWS.rglob("*.py")):
        for cls in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(cls, ast.ClassDef):
                for node in cls.body:
                    yield path, cls, node


def _location(path: Path, node: ast.stmt) -> str:
    return f"{path.relative_to(REPO)}:{node.lineno}"


def test_no_get_args():
    found = [
        _location(path, node)
        for path, _, node in _class_members()
        if isinstance(node, ast.FunctionDef) and re.match(r"_?get_args", node.name)
    ]
    assert not found, f"get_args is called get_params: {found}"


def test_input_and_params_functions_take_wildcards():
    allowed = (["self", "wildcards"], ["self", "wildcards", "input"])
    found = []
    for path, cls, node in _class_members():
        if not isinstance(node, ast.FunctionDef):
            continue
        if not re.fullmatch(r"_get_(input_files|params)_\w+", node.name):
            continue
        if (cls.name, node.name) in CHECKPOINT_INPUT_FUNCTIONS:
            continue
        args = node.args
        if [a.arg for a in args.args] not in allowed or args.vararg or args.kwarg:
            found.append(f"{_location(path, node)} {cls.name}.{node.name}")
    assert not found, f"use (self, wildcards) or (self, wildcards, input): {found}"


def test_actions_are_tuples():
    found = [
        f"{_location(path, node)} {cls.name}"
        for path, cls, node in _class_members()
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "actions" for t in node.targets)
        and not isinstance(node.value, ast.Tuple)
    ]
    assert not found, f"actions must be a tuple: {found}"


def test_snakefiles_pass_input_and_params_functions():
    found = []
    rule_files = sorted([*WORKFLOWS.glob("*/Snakefile"), *WORKFLOWS.glob("*/*.rules")])
    for path in rule_files:
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            unpacked = re.search(r"\*\*\(?wf\.get_(input_files|params)\(", line)
            stray_params = "wf.get_params(" in line and "args=wf.get_params(" not in line
            if unpacked or stray_params:
                found.append(f"{path.relative_to(REPO)}:{lineno}: {line.strip()}")
    assert not found, f"use unpack(wf.get_input_files(...)) and args=wf.get_params(...): {found}"


# Path conventions ---------------------------------------------------------------------------------


def _tool_names_by_step() -> dict[str, set[str]]:
    """Return the lower-case values of each step's ``Tool`` enums."""
    result = {}
    for step in WORKFLOW_REGISTRY:
        module = importlib.import_module(f"snappy_pipeline.workflows.{step}.model")
        result[step] = {
            str(member.value).lower()
            for name, obj in vars(module).items()
            if isinstance(obj, type) and issubclass(obj, enum.Enum) and "Tool" in name
            for member in obj
        }
    return result


def _path_violations() -> dict[str, list[str]]:
    """Return step -> violations of the path conventions in all DAG snapshots.

    Paths are normalized back to templates by replacing wildcard values with ``{name}``.
    """
    tools = _tool_names_by_step()
    steps = sorted(WORKFLOW_REGISTRY, key=len, reverse=True)
    found: dict[str, set[str]] = {}
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            if job["rule"] == "all":
                continue
            step = next(s for s in steps if job["rule"].startswith(s + "_"))
            values = sorted(job["wildcards"].items(), key=lambda kv: (-len(str(kv[1])), kv[0]))
            for path in job["output"] + job["log"]:
                template = re.sub(r"^tasks/[^/]+/", "", path)
                for name, value in values:
                    template = template.replace(str(value), "{%s}" % name)
                tokens = set(re.split(r"[/._\-]", re.sub(r"\{\w+\}", "", template)))
                if tokens & tools[step]:
                    found.setdefault(step, set()).add(f"tool name: {template}")
                if re.search(r"\}_", template):
                    found.setdefault(step, set()).add(f"underscore after wildcard: {template}")
    return {step: sorted(violations) for step, violations in sorted(found.items())}


def test_path_conventions_only_improve():
    """No tool names in paths and ``{entity}.{suffix}`` naming (plans.md C3)."""
    current = _path_violations()
    if os.environ.get("SNAPPY_UPDATE_SNAPSHOTS"):
        PATH_BASELINE.write_text(json.dumps(current, indent=1) + "\n", encoding="utf-8")
        return

    baseline = json.loads(PATH_BASELINE.read_text(encoding="utf-8"))
    flat_current = {(step, v) for step, vs in current.items() for v in vs}
    flat_baseline = {(step, v) for step, vs in baseline.items() for v in vs}
    new = sorted(flat_current - flat_baseline)
    fixed = sorted(flat_baseline - flat_current)
    assert not new, f"new path convention violations: {new}"
    assert not fixed, (
        f"{len(fixed)} violations are fixed; shrink the baseline with SNAPPY_UPDATE_SNAPSHOTS=1"
    )


def test_step_snakefiles_fetch_their_workflow_object():
    found = []
    for path in sorted(WORKFLOWS.glob("*/Snakefile")):
        text = path.read_text(encoding="utf-8")
        if path.parent.name == "link_in":
            continue  # configuration carrier without rules
        if 'wf = task_instance(config["__task_name__"])' not in text or re.search(
            r"Workflow\(", text
        ):
            found.append(str(path.relative_to(REPO)))
    assert not found, f"use wf = task_instance(config['__task_name__']) only: {found}"


#: Allowed ``depends_on`` keys: named after the data, with role prefixes where a step needs two
#: inputs of one kind. ``link_in`` remains for the external-file export steps until plans.md F1.
DEPENDS_ON_KEYS = {
    "reads",
    "alignments",
    "variants",
    "somatic_variants",
    "germline_variants",
    "combined_variants",
    "annotated_variants",
    "phased_variants",
    "copy_number",
    "structural_variants",
    "fusions",
    "expression",
    "hla_types",
    "strandedness",
    "panel_of_normals",
    "reference",
    "index",
    "link_in",
}


def test_depends_on_keys_use_the_vocabulary():
    found = []
    for step, workflow_cls in WORKFLOW_REGISTRY.items():
        depends_on = workflow_cls.config_model_class.model_fields.get("depends_on")
        if depends_on is not None:
            keys = set(depends_on.annotation.model_fields) - DEPENDS_ON_KEYS
            found += [f"{step}.{key}" for key in sorted(keys)]
    assert not found, f"depends_on keys outside the vocabulary in docs/dev_conventions.rst: {found}"
