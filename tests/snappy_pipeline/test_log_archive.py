# -*- coding: utf-8 -*-
"""Tests for the log archive of a run in ``snappy_pipeline.log_archive``."""

import contextlib
import shutil
import sqlite3
import tarfile
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main
from snappy_pipeline.log_archive import (
    archive_after_run,
    build_archive,
    git_state,
    run_log_files,
)
from snappy_pipeline.orchestration import register_hooks

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")

SNAKEMAKE_LOG = """\
Building DAG of jobs...
[Sat Oct 10 12:00:00 2026]
rule mapping_bwa_run:
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
    # closing() closes the connection; the connection's own context manager only commits
    with contextlib.closing(sqlite3.connect(log_dir / "snkmt.sqlite")) as db:
        db.execute("create table jobs (id integer)")
        db.commit()
    return snakemake_log


def test_run_log_files_lists_job_and_slurm_logs(tmp_path):
    assert run_log_files(_project(tmp_path)) == [
        "tasks/mapping/work/L1/log/L1.log",
        "tasks/mapping/work/L1/log/L1.conda_info.txt",
        "tasks/calling/work/log/scatter.log",
        "slurm_log/mapping_run/4711.log",
    ]


def test_build_archive_holds_the_run_logs_and_the_snkmt_database(tmp_path):
    """The job logs with the scripts next to them, the slurm logs and the database."""
    snakemake_log = _project(
        tmp_path,
        "tasks/mapping/work/L1/log/L1.log",
        "tasks/mapping/work/L1/log/L1.sh",
        "tasks/mapping/work/L1/log/L1.R",
        "slurm_log/mapping_run/4711.log",
    )
    archive = build_archive(tmp_path, snakemake_log)

    assert archive == tmp_path / "logs" / f"{STAMP}.tar.gz"
    with tarfile.open(archive) as tar:
        assert sorted(tar.getnames()) == [
            f".snakemake/log/{STAMP}.snakemake.log",
            ".snakemake/log/snkmt.sqlite",
            "slurm_log/mapping_run/4711.log",
            "tasks/mapping/work/L1/log/L1.R",
            "tasks/mapping/work/L1/log/L1.log",
            "tasks/mapping/work/L1/log/L1.sh",
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


def _hooks(workflow_hooks: dict) -> SimpleNamespace:
    return SimpleNamespace(
        onstart=lambda func: workflow_hooks.setdefault("start", func),
        onsuccess=lambda func: workflow_hooks.setdefault("success", func),
        onerror=lambda func: workflow_hooks.setdefault("error", func),
    )


def _repository(tmp_path, git) -> str:
    """A repository with config.yaml committed; returns the commit."""
    git(tmp_path, "init", "-q")
    (tmp_path / "config.yaml").write_text("tasks: []\n")
    git(tmp_path, "add", "config.yaml")
    git(tmp_path, "commit", "-q", "-m", "Initialize the snappy project")
    return git(tmp_path, "rev-parse", "HEAD").strip()


@needs_git
def test_git_state_records_the_commit_and_the_uncommitted_changes(tmp_path, git):
    commit = _repository(tmp_path, git)
    state = git_state(tmp_path)
    assert (state.commit, state.status, state.diff) == (commit, "", "")
    assert state.text == f"commit {commit}\n"

    (tmp_path / "config.yaml").write_text("tasks: [a]\n")
    state = git_state(tmp_path)
    assert state.changed_files == ["config.yaml"]
    assert state.text.startswith(f"commit {commit}\nuncommitted changes:\n M config.yaml\n")
    assert "+tasks: [a]" in state.diff


def test_git_state_is_none_outside_a_repository(tmp_path):
    assert git_state(tmp_path) is None


def test_build_archive_holds_the_git_record(tmp_path):
    snakemake_log = _project(tmp_path)
    (tmp_path / ".snakemake" / "log" / f"{STAMP}.git.txt").write_text("commit abc\n")
    with tarfile.open(build_archive(tmp_path, snakemake_log)) as tar:
        assert f".snakemake/log/{STAMP}.git.txt" in tar.getnames()


@needs_git
def test_register_hooks_records_the_git_state_when_the_run_starts(
    tmp_path, git, monkeypatch, capsys
):
    hooks = {}
    register_hooks(_hooks(hooks))
    commit = _repository(tmp_path, git)
    (tmp_path / "config.yaml").write_text("tasks: [a]\n")
    snakemake_log = _project(tmp_path)
    monkeypatch.chdir(tmp_path)

    hooks["start"]([str(snakemake_log)])

    record = (tmp_path / ".snakemake" / "log" / f"{STAMP}.git.txt").read_text()
    assert record.startswith(f"commit {commit}\nuncommitted changes:\n")
    assert "[snappy] Warning: uncommitted changes in config.yaml" in capsys.readouterr().err


def test_register_hooks_archives_the_logs_when_the_run_ends(tmp_path, monkeypatch, capsys):
    hooks = {}
    workflow = _hooks(hooks)
    register_hooks(workflow)
    snakemake_log = _project(tmp_path)
    monkeypatch.chdir(tmp_path)

    hooks["error"]([str(snakemake_log)])

    assert (tmp_path / "logs" / f"{STAMP}.tar.gz").is_file()
    assert "Something went wrong" in capsys.readouterr().err
