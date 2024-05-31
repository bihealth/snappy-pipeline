# -*- coding: utf-8 -*-
"""Shared fixtures for the apps unit tests"""

import os.path
import textwrap

import pytest

import snappy_pipeline.apps
import snappy_pipeline.workflows
from tests.snappy_pipeline.workflows.conftest import (
    fake_fs,
    germline_sheet_fake_fs,
    germline_sheet_tsv,
)


@pytest.fixture(scope="module")  # otherwise: performance issues
def germline_sheet_config_yaml():
    """Return YAML parsing result for (germline) configuration"""
    return textwrap.dedent(
        r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

        step_config: {}

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
        """
    ).lstrip()


@pytest.fixture(scope="module")  # otherwise: performance issues
def germline_sheet_ngs_mapping_config_yaml():
    """Pipeline step ``config.yaml`` file."""
    return textwrap.dedent(
        r"""
        pipeline_step:
          name: ngs_mapping
          version: 1

        $ref: 'file://../.snappy_pipeline/config.yaml'
        """
    ).lstrip()


@pytest.fixture
def germline_sheet_fake_noproject_fs(
    fake_fs, germline_sheet_tsv, germline_sheet_fake_fs, germline_sheet_config_yaml
):
    # Create /dev/null
    if not fake_fs.os.path.exists("/dev/null"):
        fake_fs.fs.create_file("/dev/null")
    # Make the templates from apps directory visible
    fake_fs.fs.add_real_directory(os.path.dirname(snappy_pipeline.apps.__file__))
    # Create workspace (containing projects) and go there
    fake_fs.fs.create_dir("/projects")
    fake_fs.os.chdir("/projects")
    return fake_fs


@pytest.fixture
def germline_sheet_fake_project_fs(
    fake_fs, germline_sheet_tsv, germline_sheet_fake_fs, germline_sheet_config_yaml
):
    fake_fs.fs.create_dir("/project-dir/.snappy_pipeline")
    # Create the configuration YAML file
    fake_fs.fs.create_file(
        "/project-dir/.snappy_pipeline/config.yaml", contents=germline_sheet_config_yaml
    )
    # Create the sample TSV file
    fake_fs.fs.create_file(
        "/project-dir/.snappy_pipeline/sheet.tsv",
        contents=germline_sheet_tsv,
        create_missing_dirs=True,
    )
    # Create /dev/null
    if not fake_fs.os.path.exists("/dev/null"):
        fake_fs.fs.create_file("/dev/null", create_missing_dirs=True)
    # Make the templates from apps directory visible
    fake_fs.fs.add_real_directory(os.path.dirname(snappy_pipeline.apps.__file__))
    # Go into project
    fake_fs.os.chdir("/project-dir")
    return fake_fs


@pytest.fixture
def germline_sheet_fake_project_ngs_mapping_fs(
    germline_sheet_fake_project_fs, germline_sheet_ngs_mapping_config_yaml
):
    fake_fs = germline_sheet_fake_project_fs
    fake_fs.fs.create_dir("/project-dir/ngs_mapping")
    fake_fs.fs.create_file(
        "/project-dir/ngs_mapping/config.yaml",
        contents=germline_sheet_ngs_mapping_config_yaml,
        create_missing_dirs=True,
    )
    # Make the snappy_pipeline workflows visible
    fake_fs.fs.add_real_directory(snappy_pipeline.workflows.__path__._path[0])
    # Go into pipeline step
    fake_fs.os.chdir("/project-dir/ngs_mapping")
    return fake_fs
