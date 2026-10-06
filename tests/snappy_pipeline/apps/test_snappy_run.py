# -*- coding: utf-8 -*-
"""Tests for the ``snappy run`` subcommand."""

import os

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


def test_run_lets_user_arguments_override_snkmt_database(mocker, tmp_path):
    argv = _invoke_run(
        mocker, ["--directory", str(tmp_path), "--", "--logger-snkmt-db", "/custom.sqlite"]
    )

    # Snakemake uses the last occurrence of an option, so user arguments must come last.
    assert argv[-2:] == ["--logger-snkmt-db", "/custom.sqlite"]
    assert argv.index("--logger-snkmt-db") < len(argv) - 2
