# -*- coding: utf-8 -*-
"""Tests for ``snappy-snake app``"""

import os.path

import pytest

from snappy_pipeline.apps import snappy_snake
import snappy_pipeline.workflows
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_snappy_snake_help(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check whether the call to ``snappy-snake --help`` works."""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_snake", fake_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", fake_fs, mocker)
    m = mocker.MagicMock()
    mocker.patch("snappy_pipeline.apps.snappy_snake.snakemake_main", m)
    # Run the code under test
    with pytest.raises(SystemExit) as excinfo:
        snappy_snake.main(["--help", "--verbose"])
    assert "0" == str(excinfo.value)
    # Check assersions
    m.assert_not_called()


def test_snappy_snake_list_output(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check whether the call to ``snappy-snake -S`` works."""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_snake", fake_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", fake_fs, mocker)
    m = mocker.MagicMock(return_value=0)
    mocker.patch("snappy_pipeline.apps.snappy_snake.snakemake_main", m)
    # Run the code under test
    assert 0 == snappy_snake.main(["-S", "--verbose"])
    # Check assersions
    p = os.path.realpath(snappy_pipeline.workflows.__path__._path[0] + "/..")
    m.assert_called_once_with(
        [
            "--directory",
            "/project-dir/ngs_mapping",
            "--snakefile",
            p + "/workflows/ngs_mapping/Snakefile",
            "--jobscript",
            p + "/apps/tpls/jobscript.sh",
            "--rerun-triggers",
            "mtime",
            "params",
            "input",
            "-S",
            "--verbose",
            "--cores",
            "1",
            "--software-deployment-method",
            "conda",
            "--conda-frontend",
            "mamba",
        ]
    )
