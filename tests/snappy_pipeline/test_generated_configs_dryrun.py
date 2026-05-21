from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

TASK_NAMES = sorted(WORKFLOW_REGISTRY)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False)


def _tail(text: str, n: int = 40) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


@pytest.fixture(scope="session")
def generated_task_config(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    root = _repo_root()
    out_dir = tmp_path_factory.mktemp("generated-task-configs")

    gen = _run(
        [sys.executable, "scripts/generate_task_configs.py", "--out-dir", str(out_dir)], cwd=root
    )
    gen_output = (gen.stdout or "") + "\n" + (gen.stderr or "")
    if gen.returncode != 0:
        raise AssertionError(
            f"generate_task_configs.py failed\nstdout/stderr tail:\n{_tail(gen_output)}"
        )

    return {
        "root": root,
        "out_dir": out_dir,
        "config_path": out_dir / "all-workflows" / "config.yaml",
        "generate_output": gen_output,
    }


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("task_name", TASK_NAMES, ids=TASK_NAMES)
def test_generated_config_task_closure_passes(
    task_name: str, generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    root = generated_task_config["root"]
    config_path = generated_task_config["config_path"]
    out_dir = tmp_path / "dryrun"
    dry = _run(
        [
            sys.executable,
            "scripts/dryrun_generated_configs.py",
            "--config",
            str(config_path),
            "--out-dir",
            str(out_dir),
            "--task",
            task_name,
        ],
        cwd=root,
    )
    dry_output = (dry.stdout or "") + "\n" + (dry.stderr or "")
    assert dry.returncode == 0, f"task closure failed for {task_name}\n{_tail(dry_output)}"


@pytest.mark.integration
def test_generated_config_audit_regression_guard(generated_task_config: dict[str, Any]) -> None:
    gen_output = generated_task_config["generate_output"]
    unresolved = re.findall(r"unresolved validation errors, kept best-effort config", gen_output)

    cfg_path = generated_task_config["config_path"]
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
