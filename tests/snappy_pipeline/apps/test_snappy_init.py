# -*- coding: utf-8 -*-
"""Tests for ``snappy init`` subcommand."""

from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_init_with_positional_directory(germline_sheet_fake_noproject_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["init", "project"])
    assert result.exit_code == 0

    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/README.md")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/config.yaml")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/pipeline_job.sh")


def test_init_default_cwd(germline_sheet_fake_noproject_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["init"])
    assert result.exit_code == 0

    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/README.md")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/config.yaml")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/pipeline_job.sh")
