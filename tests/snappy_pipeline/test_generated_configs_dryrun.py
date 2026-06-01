from __future__ import annotations

import os
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import Any

import pytest
import yaml

from tests.scripts.generate_task_configs import build_all_tasks, load_yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _get_task_names() -> list[str]:
    root = _repo_root()
    base_config_path = root / "tests/snappy_pipeline/fixtures/base_config.yaml"
    base_config = load_yaml(base_config_path)
    # Redirect stdout to suppress print logs during test collection
    import io
    import sys as sys_orig

    f = io.StringIO()
    sys_orig.stdout = f
    try:
        tasks = build_all_tasks(base_config, base_config_path)
    finally:
        sys_orig.stdout = sys_orig.__stdout__
    return [t["name"] for t in tasks]


TASK_NAMES = _get_task_names()


RNA_TASKS = {"ngs_mapping_star", "somatic_gene_fusion_calling_arriba"}


def _fixture_dir() -> Path:
    return _repo_root() / "tests" / "snappy_pipeline" / "fixtures"


def _task_sample_sheet(task_name: str) -> Path:
    fixture_name = "samplesheet_rna.tsv" if task_name in RNA_TASKS else "samplesheet.tsv"
    return _fixture_dir() / fixture_name


def _task_raw_folders(task_name: str) -> tuple[str, ...]:
    if task_name in RNA_TASKS:
        return ("case001subregion-T1-RNA1-mRNA_seq1",)
    return ("case001subregion-N1-DNA1-WES1", "case001subregion-T1-DNA1-WES1")


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    if "PYTHONPATH" not in env:
        env["PYTHONPATH"] = str(_repo_root())
    else:
        env["PYTHONPATH"] = str(_repo_root()) + os.pathsep + env["PYTHONPATH"]
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False, env=env)


def _tail(text: str, n: int = 40) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def to_plain_obj(obj: Any) -> Any:
    """Convert ruamel/pydantic/path-like values to plain YAML-safe builtins."""
    if obj is None:
        return None
    if isinstance(obj, bool):
        return bool(obj)
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return float(obj)
    if isinstance(obj, str):
        return str(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, Enum):
        return to_plain_obj(obj.value)
    if isinstance(obj, Mapping):
        return {str(k): to_plain_obj(v) for k, v in obj.items()}
    if isinstance(obj, Sequence) and not isinstance(obj, (str, bytes, bytearray)):
        return [to_plain_obj(v) for v in obj]
    if isinstance(obj, set):
        return [to_plain_obj(v) for v in sorted(obj, key=lambda x: str(x))]
    return str(obj)


@pytest.fixture(scope="session")
def generated_task_config(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    root = _repo_root()
    out_dir = tmp_path_factory.mktemp("generated-task-configs")

    gen = _run(
        [sys.executable, "tests/scripts/generate_task_configs.py", "--out-dir", str(out_dir)],
        cwd=root,
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


def dependency_closure(task_name: str, tasks_by_name: dict[str, dict[str, Any]]) -> set[str]:
    seen = set()
    stack = [task_name]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        task = tasks_by_name.get(current)
        if not task:
            continue
        task_config = task.get("config", {}) if isinstance(task, dict) else {}
        depends_on = task_config.get("depends_on", {}) if isinstance(task_config, dict) else {}
        if not isinstance(depends_on, dict):
            continue
        for dep_task_name in depends_on.values():
            if isinstance(dep_task_name, str) and dep_task_name and dep_task_name not in seen:
                stack.append(dep_task_name)
    return seen


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("task_name", TASK_NAMES, ids=TASK_NAMES)
def test_generated_config_task_closure_passes(
    task_name: str, generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    root = generated_task_config["root"]
    config_path = generated_task_config["config_path"]

    # Load config.yaml
    config = load_yaml(config_path)
    tasks = config.get("tasks", [])
    tasks_by_name = {t["name"]: t for t in tasks if isinstance(t, dict) and "name" in t}

    # Calculate closure
    closure = dependency_closure(task_name, tasks_by_name)
    tasks_subset = [t for t in tasks if isinstance(t, dict) and t.get("name") in closure]

    # Construct closure config
    closure_config = {
        "static_data_config": config.get("static_data_config", {}),
        "tasks": tasks_subset,
        "data_sets": config.get("data_sets", {}),
    }

    # Redirect search_paths to a temporary raw directory inside the tmp_path
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    if "data_sets" in closure_config:
        for ds_name, ds_config in closure_config["data_sets"].items():
            if isinstance(ds_config, dict):
                ds_config["file"] = str(_task_sample_sheet(task_name))
                ds_config["search_paths"] = [str(raw_dir)]

    # Touch dummy FASTQ files
    for folder in _task_raw_folders(task_name):
        f_dir = raw_dir / folder
        f_dir.mkdir(parents=True, exist_ok=True)
        (f_dir / f"{folder}.R1.fastq.gz").touch()
        (f_dir / f"{folder}.R2.fastq.gz").touch()
        (raw_dir / f"{folder}.R1.fastq.gz").touch()
        (raw_dir / f"{folder}.R2.fastq.gz").touch()

    # Write config.yaml directly in tmp_path (no .snappy_pipeline subfolder!)
    closure_config_path = tmp_path / "config.yaml"
    closure_config_plain = to_plain_obj(closure_config)
    with closure_config_path.open("wt", encoding="utf-8") as f:
        yaml.safe_dump(closure_config_plain, f, sort_keys=False)

    # Ensure we emitted parseable YAML before invoking snappy/snakemake.
    try:
        reloaded = yaml.safe_load(closure_config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise AssertionError(
            f"Generated closure config is invalid YAML for {task_name}: {e}"
        ) from e
    assert isinstance(reloaded, dict), f"Generated closure config is not a mapping for {task_name}"

    # Run the dryrun command
    cmd = [
        sys.executable,
        "-m",
        "snappy_pipeline.apps.snappy_cli",
        "run",
        "-d",
        str(tmp_path),
        "--task",
        task_name,
        "--",
        "-n",
        "--cores",
        "1",
    ]
    dry = _run(cmd, cwd=root)
    dry_output = (dry.stdout or "") + "\n" + (dry.stderr or "")

    assert dry.returncode == 0, (
        f"task closure dryrun failed for {task_name}\nstdout/stderr excerpt:\n{dry_output}"
    )


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
    assert len(unresolved) <= int(os.environ.get("SNAPPY_MAX_UNRESOLVED_CONFIGS", "35"))
    assert auto_count <= int(os.environ.get("SNAPPY_MAX_AUTO_PLACEHOLDERS", "70"))
