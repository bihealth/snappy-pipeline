# -*- coding: utf-8 -*-
"""Tests for the log archive of a run in ``snappy_pipeline.log_archive``."""

import sqlite3
import tarfile
from types import SimpleNamespace

from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main
from snappy_pipeline.log_archive import archive_after_run, build_archive, run_log_files
from snappy_pipeline.orchestration import register_hooks

SNAKEMAKE_LOG = """\
Building DAG of jobs...
[Sat Oct 10 12:00:00 2026]
rule mapping_ngs_mapping_bwa_run:
    input: raw/L1.R1.fastq.gz
    output: tasks/mapping/work/L1/out/L1.bam
    log: tasks/mapping/work/L1/log/L1.log, tasks/mapping/work/L1/log/L1.conda_info.txt
    jobid: 3
Job 3 has been submitted with SLURM jobid 4711 (log: slurm_log/mapping_run/4711.log).
localrule calling_scatter:
    output: tasks/calling/work/scatter/1-of-2.region.bed
    log: tasks/calling/work/log/scatter.log
    jobid: 4
"""

STAMP = "2026-10-10T120000.000000"


def _project(tmp_path, *files):
    log_dir = tmp_path / ".snakemake" / "log"
    log_dir.mkdir(parents=True)
    snakemake_log = log_dir / f"{STAMP}.snakemake.log"
    snakemake_log.write_text(SNAKEMAKE_LOG)
    for name in files:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text(name)
    with sqlite3.connect(log_dir / "snkmt.sqlite") as db:
        db.execute("create table jobs (id integer)")
    return snakemake_log


def test_run_log_files_lists_job_and_slurm_logs(tmp_path):
    assert run_log_files(_project(tmp_path)) == [
        "tasks/mapping/work/L1/log/L1.log",
        "tasks/mapping/work/L1/log/L1.conda_info.txt",
        "tasks/calling/work/log/scatter.log",
        "slurm_log/mapping_run/4711.log",
    ]


def test_build_archive_holds_the_run_logs_and_the_snkmt_database(tmp_path):
    snakemake_log = _project(
        tmp_path, "tasks/mapping/work/L1/log/L1.log", "slurm_log/mapping_run/4711.log"
    )
    archive = build_archive(tmp_path, snakemake_log)

    assert archive == tmp_path / "logs" / f"{STAMP}.tar.gz"
    with tarfile.open(archive) as tar:
        assert sorted(tar.getnames()) == [
            f".snakemake/log/{STAMP}.snakemake.log",
            ".snakemake/log/snkmt.sqlite",
            "slurm_log/mapping_run/4711.log",
            "tasks/mapping/work/L1/log/L1.log",
        ]


def test_archive_after_run_reports_failures_without_raising(tmp_path, capsys):
    archive_after_run(tmp_path, [str(tmp_path / "missing.snakemake.log")])
    assert "Could not archive the logs of this run" in capsys.readouterr().err


def test_snappy_logs_archives_the_newest_run(tmp_path):
    _project(tmp_path, "tasks/calling/work/log/scatter.log")
    result = CliRunner().invoke(main, ["logs", "--directory", str(tmp_path)])

    assert result.exit_code == 0, result.output
    with tarfile.open(tmp_path / "logs" / f"{STAMP}.tar.gz") as tar:
        assert "tasks/calling/work/log/scatter.log" in tar.getnames()


def test_register_hooks_archives_the_logs_when_the_run_ends(tmp_path, monkeypatch, capsys):
    hooks = {}
    workflow = SimpleNamespace(
        onsuccess=lambda func: hooks.setdefault("success", func),
        onerror=lambda func: hooks.setdefault("error", func),
    )
    register_hooks(workflow)
    snakemake_log = _project(tmp_path)
    monkeypatch.chdir(tmp_path)

    hooks["error"]([str(snakemake_log)])

    assert (tmp_path / "logs" / f"{STAMP}.tar.gz").is_file()
    assert "Something went wrong" in capsys.readouterr().err
