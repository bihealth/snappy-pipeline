# -*- coding: utf-8 -*-
"""Tests for ``snappy-start-step`` app."""

from snappy_pipeline.apps import snappy_start_step

from ..workflows.conftest import patch_module_fs


def test_start_step_ngs_mapping(germline_sheet_fake_project_fs, mocker):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.apps.snappy_start_step", germline_sheet_fake_project_fs, mocker
    )
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_project_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_project_fs, mocker)
    # Run the code under test
    assert snappy_start_step.main(["--step", "ngs_mapping"]) is None  # and no exception
    # Check result
    germline_sheet_fake_project_fs.os.path.exists("/project-dir/ngs_mapping/config.yaml")
    germline_sheet_fake_project_fs.os.path.exists("/project-dir/ngs_mapping/pipeline_job.sh")
    germline_sheet_fake_project_fs.os.path.exists(
        "/project-dir/.snappy_pipeline/config.ruamel_yaml.bak"
    )
