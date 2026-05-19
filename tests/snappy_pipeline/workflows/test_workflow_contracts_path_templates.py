# -*- coding: utf-8 -*-
"""Central policy checks for workflow path template naming.

These checks focus on workflows that were already migrated to external task multiplexing.
They ensure path templates do not encode config-derived tool/caller dimensions.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

WORKFLOW_ROOT = Path("snappy_pipeline/workflows")

ALL_POLICY_STEPS = sorted(
    step_name
    for step_name in WORKFLOW_REGISTRY
    if (WORKFLOW_ROOT / step_name / "__init__.py").exists()
)

# Backlog snapshot: workflows still encoding config-derived naming dimensions.
# Keep this list explicit so newly introduced regressions fail immediately.
NEEDS_ADAPTATION = set()

POLICY_ENFORCED_STEPS = sorted(set(ALL_POLICY_STEPS) - NEEDS_ADAPTATION)

# Config-derived wildcard dimensions that should not appear in path templates.
FORBIDDEN_PLACEHOLDER_RE = re.compile(
    r"\{+\s*(?:"
    r"mapper|caller|var_caller|anno_caller|tool|trimmer|"
    r"mapping_tool|expression_tool|copy_number_tool"
    r")\s*\}+"
)

# Common legacy hardcoded tool/caller prefixes before library wildcards.
FORBIDDEN_PREFIX_RE = re.compile(
    r"\b(?:"
    r"bwa|bwa_mem2|minimap2|star|mbcs|bbduk|fastp|"
    r"mutect2|vep|cnvkit|sequenza|purecn"
    r")\.\{(?:library_name|tumor_library|normal_library|index_ngs_library)\}"
)


def _scan_file(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    violations: list[str] = []
    in_doc_block = False

    for i, line in enumerate(text.splitlines(), start=1):
        triple_quote_count = line.count('"""') + line.count("'''")
        if in_doc_block:
            if triple_quote_count % 2 == 1:
                in_doc_block = False
            continue
        if triple_quote_count % 2 == 1:
            in_doc_block = True
            continue

        # Ignore comments/doc lines to reduce noise; enforce actual code templates.
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if FORBIDDEN_PLACEHOLDER_RE.search(line):
            violations.append(f"{path}:{i}: forbidden wildcard in template: {stripped}")
        if FORBIDDEN_PREFIX_RE.search(line):
            violations.append(f"{path}:{i}: forbidden tool/caller prefix in template: {stripped}")

    return violations


@pytest.mark.parametrize("step_name", POLICY_ENFORCED_STEPS)
def test_no_config_derived_dimensions_in_path_templates(step_name: str):
    init_py = WORKFLOW_ROOT / step_name / "__init__.py"
    assert init_py.exists(), f"Missing workflow implementation file: {init_py}"

    violations = _scan_file(init_py)
    assert not violations, "\n".join(violations)


def test_path_template_policy_backlog_snapshot():
    offenders = set()
    for step_name in ALL_POLICY_STEPS:
        init_py = WORKFLOW_ROOT / step_name / "__init__.py"
        if _scan_file(init_py):
            offenders.add(step_name)

    assert offenders == NEEDS_ADAPTATION
