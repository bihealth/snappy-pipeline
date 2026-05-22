# -*- coding: utf-8 -*-
"""Tests for ``snappy run`` subcommand."""

import os.path
from click.testing import CliRunner
import pytest

import snappy_pipeline.workflows
from snappy_pipeline.apps import snappy_cli
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_snappy_snake_help(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check whether the call to ``snappy run --help`` works."""
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_cli", fake_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", fake_fs, mocker)
    m = mocker.MagicMock()
    mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", m)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["run", "--help"])
    assert result.exit_code == 0
    m.assert_not_called()


def test_snappy_snake_list_output(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check whether the call to ``snappy run`` works."""
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_cli", fake_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", fake_fs, mocker)
    m = mocker.MagicMock(return_value=0)
    mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", m)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["run", "--verbose"])
    assert result.exit_code == 0

    p = os.path.realpath(snappy_pipeline.workflows.__path__[0] + "/..")
    m.assert_called_once_with(
        [
            "--directory",
            "/project-dir/ngs_mapping",
            "--snakefile",
            p + "/Snakefile",
            "--config",
            "dump_orchestrator=True",
        ]
    )
