#!/usr/bin/env python3
"""Build the orchestrator DAG for one task and dump all its jobs as JSON.

For every job, the dump records the rule, wildcards, input, output and log files, evaluated
params, threads, resources and the wrapper. It is the basis of the DAG snapshot tests in
``tests/snappy_pipeline/test_generated_configs_dryrun.py``. It fails if a rule with a wrapper
has no log named ``log``.

Usage: dump_dag.py --directory <project dir with config.yaml> --task <task name> --output <json>
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from snakemake.api import SnakemakeApi
from snakemake.io import is_flagged
from snakemake.settings.types import ConfigSettings, DAGSettings, OutputSettings, ResourceSettings

import snappy_pipeline

#: Resources that depend on the machine running the test, not on the workflow.
_MACHINE_RESOURCES = {"tmpdir", "shell_exec"}


def _plain(value: Any) -> Any:
    """Convert Snakemake and Python containers into JSON-serializable builtins."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return sorted((_plain(v) for v in value), key=str)
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return str(value)


def _job_record(job) -> dict[str, Any]:
    record = {
        "rule": job.rule.name,
        "wildcards": _plain(dict(job.wildcards_dict)),
        "input": [str(f) for f in job.input],
        "output": [str(f) for f in job.output],
        "log": [str(f) for f in job.log],
        "benchmark": str(job.benchmark) if job.benchmark else None,
        "temp": sorted(str(f) for f in job.temp_output),
        "params": _plain(dict(job.params.items())),
        "threads": job.threads,
        "resources": {
            k: _plain(v)
            for k, v in job.resources.items()
            if not k.startswith("_") and k not in _MACHINE_RESOURCES
        },
        "wrapper": job.rule.wrapper,
    }
    # Outputs that the Snakemake report shows; only when there are some, as most jobs have none
    if report := sorted(str(f) for f in job.output if is_flagged(f, "report")):
        record["report"] = report
    # directory() outputs, which have no md5 file; only when there are some, like the report
    if directories := sorted(str(f) for f in job.output if is_flagged(f, "directory")):
        record["directory"] = directories
    # The image of rules that run in a container (``software: container(...)``)
    spec = getattr(job.rule.software_env_specs, "software_spec", None)
    if getattr(spec, "kind", None) == "container":
        record["container"] = spec.image_uri
    # Rules that run on the submit node instead of the cluster (``localrule: True``)
    if job.rule.name in job.rule.workflow._localrules:
        record["local"] = True
    # The job group as configured: Rule.group is None outside cluster runs and dry runs
    if job.rule._group is not None:
        record["group"] = job.rule._group
    return record


def dump_dag(directory: Path, task: str, frozen: bool = False) -> list[dict[str, Any]]:
    """Return the job records of the DAG that ``snappy run --task <task> [--frozen]`` would build."""
    snakefile = Path(snappy_pipeline.__file__).parent / "Snakefile"
    with SnakemakeApi(OutputSettings()) as api:
        workflow_api = api.workflow(
            resource_settings=ResourceSettings(cores=1),
            config_settings=ConfigSettings(config={"task": task, "frozen": frozen}),
            snakefile=snakefile,
            workdir=directory,
        )
        dag_api = workflow_api.dag(dag_settings=DAGSettings(forceall=True))
        # Snakemake has no public API that returns the DAG's jobs; this mirrors printdag().
        workflow = dag_api.workflow_api._workflow
        workflow._prepare_dag(forceall=True, ignore_incomplete=True, lock_warn_only=True)
        workflow._build_dag()
        records = [_job_record(job) for job in workflow.dag.jobs]
        # The wrapper base classes fail at run time without snakemake.log.log.
        unnamed = sorted(
            {
                job.rule.name
                for job in workflow.dag.jobs
                if job.rule.wrapper and "log" not in job.log.keys()
            }
        )
        if unnamed:
            raise SystemExit(f"wrapper rules without a log named 'log': {', '.join(unnamed)}")
    return sorted(records, key=lambda r: (r["rule"], json.dumps(r["wildcards"]), r["output"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--task", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frozen", action="store_true")
    args = parser.parse_args()

    records = dump_dag(args.directory.resolve(), args.task, args.frozen)
    args.output.write_text(json.dumps(records, indent=1, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
