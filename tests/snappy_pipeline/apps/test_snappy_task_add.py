# -*- coding: utf-8 -*-
"""Tests for ``snappy task add`` subcommand."""

from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli
from ..workflows.conftest import patch_module_fs


def test_task_add_ngs_mapping(germline_sheet_fake_project_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "add", "ngs_mapping"])
    assert result.exit_code == 0

    assert germline_sheet_fake_project_fs.os.path.exists("/project-dir/pipeline_job.sh")
    assert germline_sheet_fake_project_fs.os.path.exists("/project-dir/config.yaml.bak")
    with germline_sheet_fake_project_fs.open("/project-dir/config.yaml", "rt") as f:
        content = f.read()
        assert "ngs_mapping" in content


def test_task_add_custom_task_name(germline_sheet_fake_project_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "add", "custom_mapping=ngs_mapping"])
    assert result.exit_code == 0

    assert germline_sheet_fake_project_fs.os.path.exists("/project-dir/pipeline_job.sh")
    assert germline_sheet_fake_project_fs.os.path.exists("/project-dir/config.yaml.bak")
    with germline_sheet_fake_project_fs.open("/project-dir/config.yaml", "rt") as f:
        content = f.read()
        assert "custom_mapping" in content
        assert "step: ngs_mapping" in content


def test_task_add_config_mode_full(germline_sheet_fake_project_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "add", "ngs_mapping", "--config-mode", "full"])
    assert result.exit_code == 0
    with germline_sheet_fake_project_fs.open("/project-dir/config.yaml", "rt") as f:
        content = f.read()
        assert "tool: bwa" in content or "tool:" in content


def test_task_add_config_mode_minimal(germline_sheet_fake_project_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(
        snappy_cli.main, ["task", "add", "ngs_mapping", "--config-mode", "minimal"]
    )
    assert result.exit_code == 0
    with germline_sheet_fake_project_fs.open("/project-dir/config.yaml", "rt") as f:
        content = f.read()
        assert "step: ngs_mapping" in content
