from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False)


def _tail(text: str, n: int = 40) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def _run_generate_and_dryrun() -> tuple[Path, str]:
    root = _repo_root()

    gen = _run([sys.executable, "scripts/generate_task_configs.py"], cwd=root)
    if gen.returncode != 0:
        raise AssertionError(
            "generate_task_configs.py failed\n"
            f"stdout/stderr tail:\n{_tail((gen.stdout or '') + '\n' + (gen.stderr or ''))}"
        )

    dry = _run([sys.executable, "scripts/dryrun_generated_configs.py"], cwd=root)
    dry_output = (dry.stdout or "") + "\n" + (dry.stderr or "")
    if dry.returncode != 0:
        raise AssertionError(
            f"dryrun_generated_configs.py failed\nstdout/stderr tail:\n{_tail(dry_output)}"
        )

    report_json = root / "scratch/generated-configs/reports/dryrun_shards_report.json"
    if not report_json.exists():
        raise AssertionError(f"Expected report not found: {report_json}")

    return report_json, dry_output


@pytest.mark.integration
@pytest.mark.slow
def test_generated_config_shards_all_pass() -> None:
    report_json, _ = _run_generate_and_dryrun()
    report = json.loads(report_json.read_text(encoding="utf-8"))

    summary = report.get("summary", {})
    total = int(summary.get("total", -1))
    failed = int(summary.get("failed", -1))

    assert total >= 40
    assert failed == 0


@pytest.mark.integration
@pytest.mark.slow
def test_generated_config_audit_regression_guard() -> None:
    root = _repo_root()

    gen = _run([sys.executable, "scripts/generate_task_configs.py"], cwd=root)
    gen_output = (gen.stdout or "") + "\n" + (gen.stderr or "")
    if gen.returncode != 0:
        raise AssertionError(
            f"generate_task_configs.py failed\nstdout/stderr tail:\n{_tail(gen_output)}"
        )

    unresolved = re.findall(r"unresolved validation errors, kept best-effort config", gen_output)

    cfg_path = root / "scratch/generated-configs/all-workflows/config.yaml"
    cfg: dict[str, Any] = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}

    auto_count = 0

    def walk(v: Any) -> None:
        nonlocal auto_count
        if isinstance(v, dict):
            for vv in v.values():
                walk(vv)
        elif isinstance(v, list):
            for vv in v:
                walk(vv)
        elif isinstance(v, str) and v == "AUTO":
            auto_count += 1

    walk(cfg)

    # Guard against silent regressions: keep generated-config quality from getting worse.
    assert len(unresolved) <= int(os.environ.get("SNAPPY_MAX_UNRESOLVED_CONFIGS", "13"))
    assert auto_count <= int(os.environ.get("SNAPPY_MAX_AUTO_PLACEHOLDERS", "22"))
