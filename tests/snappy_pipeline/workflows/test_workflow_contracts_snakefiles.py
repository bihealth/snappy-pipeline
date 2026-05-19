# -*- coding: utf-8 -*-
"""Contract tests for Snakefile tool gating patterns."""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
WORKFLOWS_DIR = REPO_ROOT / "snappy_pipeline" / "workflows"


def _snakefiles():
    return sorted(path for path in WORKFLOWS_DIR.glob("*/Snakefile") if path.is_file())


def test_snakefiles_do_not_use_legacy_multi_tool_membership_checks():
    # A workflow task has exactly one selected tool, so rule gating must branch on `tool == ...`.
    forbidden = re.compile(r"^\s*if\s+[\"'][-_a-zA-Z0-9]+[\"']\s+in\s+tools\s*:")
    offenders = []
    for snakefile in _snakefiles():
        for lineno, line in enumerate(snakefile.read_text().splitlines(), start=1):
            if forbidden.search(line):
                offenders.append(f"{snakefile}:{lineno}:{line.strip()}")

    assert not offenders, "Legacy multi-tool gating found:\n" + "\n".join(offenders)


def test_snakefiles_do_not_wrap_single_tool_in_tools_list():
    forbidden = re.compile(r"^\s*tools\s*=\s*\[\s*str\(wf\.config\.tool\)\s*\]\s*$")
    offenders = []
    for snakefile in _snakefiles():
        for lineno, line in enumerate(snakefile.read_text().splitlines(), start=1):
            if forbidden.search(line):
                offenders.append(f"{snakefile}:{lineno}:{line.strip()}")

    assert not offenders, "Single-tool list wrappers found:\n" + "\n".join(offenders)
