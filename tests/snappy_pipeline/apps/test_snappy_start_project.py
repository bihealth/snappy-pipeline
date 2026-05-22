# -*- coding: utf-8 -*-
"""Tests for ``snappy init`` subcommand."""

from click.testing import CliRunner
from snappy_pipeline.apps import snappy_cli
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_start_project(germline_sheet_fake_noproject_fs, mocker):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.apps.snappy_cli", germline_sheet_fake_noproject_fs, mocker
    )
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)
    
    # Run the click CLI command under test
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["init", "--directory", "project"])
    assert result.exit_code == 0
    
    # Check result
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/README.md")
    assert germline_sheet_fake_noproject_fs.os.path.exists(
        "/projects/project/.snappy_pipeline/config.yaml"
    )
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/pipeline_job.sh")
