# -*- coding: utf-8 -*-
"""Tests for ``snappy`` CLI subcommands, completion, and task listing."""

from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_steps_match_workflow_registry():
    """STEPS used by the CLI must match the WORKFLOW_REGISTRY keys."""
    assert set(snappy_cli.STEPS) == set(WORKFLOW_REGISTRY)


def test_snappy_run_help(germline_sheet_fake_project_ngs_mapping_fs, mocker):
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


def test_snappy_completion():
    """Check whether ``snappy completion`` generates autocompletion code."""
    runner = CliRunner()
    result_zsh = runner.invoke(snappy_cli.main, ["completion", "zsh"])
    assert result_zsh.exit_code == 0
    assert "_SNAPPY_COMPLETE" in result_zsh.output

    result_bash = runner.invoke(snappy_cli.main, ["completion", "bash"])
    assert result_bash.exit_code == 0
    assert "_SNAPPY_COMPLETE" in result_bash.output


def test_snappy_task_steps():
    """Check whether ``snappy task steps`` lists available workflow steps."""
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "steps"])
    assert result.exit_code == 0
    assert "Available Workflow Steps in SNAPPY" in result.output
    assert "ngs_mapping" in result.output


def test_snappy_task_list_available():
    """Check whether ``snappy task list --available`` lists available steps."""
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "list", "--available"])
    assert result.exit_code == 0
    assert "Available Workflow Steps in SNAPPY" in result.output
    assert "ngs_mapping" in result.output


def test_task_param_type_shell_complete():
    """Check shell completion helper for TaskParamType."""
    param_type = snappy_cli.TaskParamType()
    completions = param_type.shell_complete(None, None, "ngs_")
    names = [c.value for c in completions]
    assert "ngs_mapping" in names

    completions_prefix = param_type.shell_complete(None, None, "my_custom_task=ngs_")
    names_prefix = [c.value for c in completions_prefix]
    assert "my_custom_task=ngs_mapping" in names_prefix


def test_complete_task_names(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check task name completion from project config.yaml."""
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_cli", fake_fs, mocker)

    ctx = mocker.MagicMock()
    ctx.params = {"directory": "/project-dir/ngs_mapping"}

    completions = snappy_cli._complete_task_names(ctx, None, "")
    names = [c.value for c in completions]
    assert "ngs_mapping" in names
