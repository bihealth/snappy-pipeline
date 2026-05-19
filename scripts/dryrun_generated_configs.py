#!/usr/bin/env python3
"""Create dependency-closure task shards from a generated config and dry-run each shard.

This is intended to make migration/fix work tractable: each step is tested in a
minimal config containing itself plus transitive `depends_on` tasks.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import ruamel.yaml as ruamel_yaml

yaml = ruamel_yaml.YAML()
yaml.default_flow_style = False
yaml.indent(sequence=4, offset=2)


@dataclass
class DryRunResult:
    task_name: str
    shard_dir: Path
    return_code: int
    ok: bool
    first_error_line: str
    error_excerpt: str


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


def dependency_closure(task_name: str, tasks_by_name: dict[str, dict[str, Any]]) -> set[str]:
    seen: set[str] = set()
    stack = [task_name]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        task = tasks_by_name.get(current)
        if not task:
            continue
        depends_on = task.get("depends_on", {})
        if not isinstance(depends_on, dict):
            continue
        for dep_task_name in depends_on.values():
            if isinstance(dep_task_name, str) and dep_task_name and dep_task_name not in seen:
                stack.append(dep_task_name)
    return seen


def build_shard_config(
    base_config: dict[str, Any], tasks_subset: list[dict[str, Any]]
) -> dict[str, Any]:
    return {
        "static_data_config": base_config.get("static_data_config", {}),
        "tasks": tasks_subset,
        "data_sets": base_config.get("data_sets", {}),
    }


def run_dry_run(config_dir: Path) -> tuple[int, str]:
    cmd = [
        "snappy-snake",
        "-d",
        str(config_dir),
        "--",
        "-n",
        "--cores",
        "1",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    output = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, output


def first_error_line(output: str) -> str:
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        if "ValidationError in file" in line:
            return line
        if "AttributeError in file" in line:
            return line
        if "ValueError in file" in line:
            return line
        if "TypeError in file" in line:
            return line
        if line.startswith("Error in rule"):
            return line
    for line in output.splitlines():
        line = line.strip()
        if "ERROR" in line:
            return line
    return "(no error line detected)"


def error_excerpt(output: str, max_lines: int = 8) -> str:
    lines = output.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if "ValidationError in file" in line or "AttributeError in file" in line:
            start = idx
            break
        if "ValueError in file" in line or "TypeError in file" in line:
            start = idx
            break
    if start is None:
        start = 0
    excerpt = lines[start : start + max_lines]
    return "\n".join(excerpt).strip()


def write_reports(out_dir: Path, results: list[DryRunResult]) -> None:
    reports_dir = out_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    json_path = reports_dir / "dryrun_shards_report.json"
    payload = {
        "summary": {
            "total": len(results),
            "ok": sum(1 for r in results if r.ok),
            "failed": sum(1 for r in results if not r.ok),
        },
        "results": [
            {
                "task_name": r.task_name,
                "shard_dir": str(r.shard_dir),
                "return_code": r.return_code,
                "ok": r.ok,
                "first_error_line": r.first_error_line,
                "error_excerpt": r.error_excerpt,
            }
            for r in results
        ],
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_path = reports_dir / "dryrun_shards_report.md"
    lines = []
    lines.append("# Dry-Run Shard Report")
    lines.append("")
    lines.append(f"- Total: {payload['summary']['total']}")
    lines.append(f"- OK: {payload['summary']['ok']}")
    lines.append(f"- Failed: {payload['summary']['failed']}")
    lines.append("")
    lines.append("## Failed")
    lines.append("")
    for r in results:
        if r.ok:
            continue
        lines.append(f"- `{r.task_name}`: {r.first_error_line}")
        if r.error_excerpt:
            lines.append("\n```text")
            lines.append(r.error_excerpt)
            lines.append("```")
    lines.append("")
    lines.append("## Passed")
    lines.append("")
    passed = [r.task_name for r in results if r.ok]
    if passed:
        lines.append("- " + ", ".join(f"`{name}`" for name in passed))
    else:
        lines.append("- (none)")
    lines.append("")

    md_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("scratch/generated-configs/all-workflows/config.yaml"),
        help="Generated all-workflows config to shard.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("scratch/generated-configs"),
        help="Base output dir for shards and reports.",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=0,
        help="Optional cap for number of tasks to process (0 means all).",
    )
    args = parser.parse_args()

    config_path = args.config.resolve()
    out_dir = args.out_dir.resolve()

    config = load_yaml(config_path)
    tasks = config.get("tasks", [])
    if not isinstance(tasks, list) or not tasks:
        raise ValueError("Config has no tasks to shard")

    tasks_by_name = {t.get("name"): t for t in tasks if isinstance(t, dict) and t.get("name")}

    results: list[DryRunResult] = []
    shards_root = out_dir / "shards"
    shards_root.mkdir(parents=True, exist_ok=True)

    ordered_task_names = [t["name"] for t in tasks if isinstance(t, dict) and "name" in t]
    if args.max_steps > 0:
        ordered_task_names = ordered_task_names[: args.max_steps]

    for task_name in ordered_task_names:
        closure = dependency_closure(task_name, tasks_by_name)
        tasks_subset = [t for t in tasks if isinstance(t, dict) and t.get("name") in closure]

        shard_dir = shards_root / task_name
        shard_config_path = shard_dir / "config.yaml"
        shard_config = build_shard_config(config, tasks_subset)
        dump_yaml(shard_config_path, shard_config)

        rc, output = run_dry_run(shard_dir)
        (shard_dir / "dryrun.log").write_text(output, encoding="utf-8")
        result = DryRunResult(
            task_name=task_name,
            shard_dir=shard_dir,
            return_code=rc,
            ok=(rc == 0),
            first_error_line=first_error_line(output),
            error_excerpt=error_excerpt(output),
        )
        results.append(result)
        status = "OK" if result.ok else "FAIL"
        print(f"[{status}] {task_name}: {result.first_error_line}")

    write_reports(out_dir, results)

    total = len(results)
    ok = sum(1 for r in results if r.ok)
    failed = total - ok
    print(f"Wrote reports to {out_dir / 'reports'}")
    print(f"Summary: total={total} ok={ok} failed={failed}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
