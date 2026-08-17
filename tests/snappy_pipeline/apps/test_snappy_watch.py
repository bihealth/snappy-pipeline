# -*- coding: utf-8 -*-
"""Tests for the ``snappy watch`` subcommand."""

import subprocess

from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli


def test_snkmt_console_cmd_prefers_console_script(mocker, tmp_path):
    fake_bin = str(tmp_path / "snkmt")
    mocker.patch("snappy_pipeline.apps.snappy_cli.shutil.which", return_value=fake_bin)

    cmd = snappy_cli._snkmt_console_cmd("/project/.snakemake/log/snkmt.sqlite")

    assert cmd == [
        fake_bin,
        "console",
        "--db-path",
        "/project/.snakemake/log/snkmt.sqlite",
    ]


def test_snkmt_console_cmd_falls_back_to_main(mocker):
    # No `snkmt` console script on PATH -> run snkmt.cli:main() with the
    # current interpreter.  `python -m snkmt.cli` is deliberately NOT used:
    # the module has no `__main__` guard and would silently do nothing.
    mocker.patch("snappy_pipeline.apps.snappy_cli.shutil.which", return_value=None)
    mocker.patch("snappy_pipeline.apps.snappy_cli.sys.executable", "/venv/bin/python")

    cmd = snappy_cli._snkmt_console_cmd("/project/snkmt.sqlite")

    assert cmd == [
        "/venv/bin/python",
        "-c",
        "from snkmt.cli import main; main()",
        "console",
        "--db-path",
        "/project/snkmt.sqlite",
    ]


def test_watch_invokes_snkmt_with_database(mocker, tmp_path):
    db_path = tmp_path / ".snakemake" / "log" / "snkmt.sqlite"
    db_path.parent.mkdir(parents=True)
    db_path.touch()
    fake_bin = str(tmp_path / "bin" / "snkmt")

    mocker.patch("snappy_pipeline.apps.snappy_cli.shutil.which", return_value=fake_bin)
    recorded = {}

    def fake_run(cmd, **kwargs):
        recorded["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0)

    mocker.patch("snappy_pipeline.apps.snappy_cli.subprocess.run", side_effect=fake_run)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["watch", "--directory", str(tmp_path)])

    assert result.exit_code == 0, result.output
    assert recorded["cmd"] == [fake_bin, "console", "--db-path", str(db_path)]


def test_watch_missing_database_errors(mocker, tmp_path):
    log_mock = mocker.patch("snappy_pipeline.apps.snappy_cli.log")
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["watch", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    error_calls = [c for c in log_mock.call_args_list if "snkmt database not found" in c.args[0]]
    assert len(error_calls) == 1
