# -*- coding: utf-8 -*-
"""Tests for the ``snappy run`` subcommand."""

import os

import click
from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli


def _invoke_run(mocker, args):
    snakemake_main = mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", return_value=0)
    result = CliRunner().invoke(snappy_cli.main, ["run", *args])
    assert result.exit_code == 0, result.output
    return snakemake_main.call_args.args[0]


def test_run_writes_snkmt_database_into_project(mocker, tmp_path):
    argv = _invoke_run(mocker, ["--directory", str(tmp_path)])

    db_arg = argv.index("--logger-snkmt-db")
    assert argv[db_arg + 1] == str(tmp_path / ".snakemake" / "log" / "snkmt.sqlite")


def test_run_makes_snkmt_database_path_absolute(mocker, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    argv = _invoke_run(mocker, ["--directory", "project"])

    db_arg = argv.index("--logger-snkmt-db")
    assert argv[db_arg + 1] == os.path.join(
        str(tmp_path), "project", ".snakemake", "log", "snkmt.sqlite"
    )


def test_run_passes_the_config_file_as_an_absolute_path(mocker, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    argv = _invoke_run(mocker, ["--directory", "project"])

    # load_project resolves relative config paths against the directory of this file
    assert argv[argv.index("--configfile") + 1] == str(tmp_path / "project" / "config.yaml")


def test_run_lets_user_arguments_override_snkmt_database(mocker, tmp_path):
    argv = _invoke_run(
        mocker, ["--directory", str(tmp_path), "--", "--logger-snkmt-db", "/custom.sqlite"]
    )

    # Snakemake uses the last occurrence of an option, so user arguments must come last.
    assert argv[-2:] == ["--logger-snkmt-db", "/custom.sqlite"]
    assert argv.index("--logger-snkmt-db") < len(argv) - 2


def test_run_uses_orchestrator_snakefile_and_default_profile(mocker, tmp_path):
    argv = _invoke_run(mocker, ["--directory", str(tmp_path)])

    assert argv[argv.index("--snakefile") + 1].endswith(
        os.path.join("snappy_pipeline", "Snakefile")
    )
    profiles = [argv[i + 1] for i, arg in enumerate(argv) if arg == "--workflow-profile"]
    assert [os.path.basename(p) for p in profiles] == ["profile"]
    assert "--config" not in argv


def test_run_passes_target_selection_as_config(mocker, tmp_path):
    argv = _invoke_run(mocker, ["--directory", str(tmp_path), "--task", "calling"])
    assert argv[argv.index("--config") + 1 :][:1] == ["task=calling"]

    argv = _invoke_run(mocker, ["--directory", str(tmp_path)])
    assert "--config" not in argv or not argv[argv.index("--config") + 1].startswith("task=")


def test_run_passes_frozen_mode_as_config(mocker, tmp_path):
    argv = _invoke_run(mocker, ["--directory", str(tmp_path), "--task", "annotation", "--frozen"])
    assert argv[argv.index("--config") + 1 :][:2] == ["task=annotation", "frozen=True"]


def test_run_frozen_requires_a_task(mocker, tmp_path):
    mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", return_value=0)
    result = CliRunner().invoke(snappy_cli.main, ["run", "--directory", str(tmp_path), "--frozen"])
    assert result.exit_code != 0
    # rich-click colors the options in its error box on GitHub Actions (GITHUB_ACTIONS)
    assert "--frozen requires --task" in click.unstyle(result.output)


def test_run_leaves_the_executor_to_the_profile_of_the_user(mocker, tmp_path):
    # Only the workflow profile; a profile of the system or the user passes through
    argv = _invoke_run(mocker, ["--directory", str(tmp_path), "--", "--profile", "slurm"])

    profiles = [argv[i + 1] for i, arg in enumerate(argv) if arg == "--workflow-profile"]
    assert [os.path.basename(p) for p in profiles] == ["profile"]
    assert argv[argv.index("--profile") + 1] == "slurm"
