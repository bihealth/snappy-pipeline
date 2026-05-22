# -*- coding: utf-8 -*-
"""Tests for ``snappy task add`` subcommand."""

from click.testing import CliRunner
from snappy_pipeline.apps import snappy_cli
from ..workflows.conftest import patch_module_fs


def test_start_step_ngs_mapping(germline_sheet_fake_project_fs, mocker):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.apps.snappy_cli", germline_sheet_fake_project_fs, mocker
    )
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)
    
    # Run the code under test
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "add", "--step", "ngs_mapping"])
    assert result.exit_code == 0
    
    # Check result
    assert germline_sheet_fake_project_fs.os.path.exists("/project-dir/pipeline_job.sh")
    assert germline_sheet_fake_project_fs.os.path.exists(
        "/project-dir/.snappy_pipeline/config.yaml.bak"
    )
    with germline_sheet_fake_project_fs.open("/project-dir/.snappy_pipeline/config.yaml", "rt") as f:
        content = f.read()
        assert "ngs_mapping" in content
