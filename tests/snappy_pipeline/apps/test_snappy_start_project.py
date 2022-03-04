# -*- coding: utf-8 -*-
"""Tests for ``snappy-start-step`` app."""

from snappy_pipeline.apps import snappy_start_project

from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_start_project(germline_sheet_fake_noproject_fs, mocker):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.apps.snappy_start_step", germline_sheet_fake_noproject_fs, mocker
    )
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)
    # Run the code under test
    assert snappy_start_project.main(["--directory", "project"]) is None  # and no exception
    # Check result
    germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/README.md")
    germline_sheet_fake_noproject_fs.os.path.exists(
        "/projects/project/.snappy_pipeline/config.yaml"
    )
