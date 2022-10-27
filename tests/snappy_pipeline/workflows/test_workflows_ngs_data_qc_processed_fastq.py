# -*- coding: utf-8 -*-
"""Tests for the ngs_data_qc workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.ngs_data_qc import NgsDataQcWorkflow

from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

        step_config:
          ngs_data_qc:
            path_link_in: "/preprocess"
            tools: ['fastqc']

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
            pedigree_field: pedigree_field
        """
        ).lstrip()
    )


@pytest.fixture
def ngs_data_qc(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs_path_link_in,
    aligner_indices_fake_fs,
    mocker,
):
    """Return NgsDataQcWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs_path_link_in, mocker
    )
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_data_qc", aligner_indices_fake_fs, mocker)
    # Construct the workflow object
    return NgsDataQcWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FastQcReportStepPart -------------------------------------------------------------------


def test_fastqc_step_part_get_args(ngs_data_qc):
    """Tests FastQcReportStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "num_threads": 1,
        "more_reads": [
            "work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz",
            "work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz",
        ],
    }
    # Get actual and assert
    actual = ngs_data_qc.get_args("fastqc", "run")(wildcards)
    assert actual == expected


def test_fastqc_step_part_get_input_files(ngs_data_qc):
    """Tests FastQcReportStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_data_qc.get_input_files("fastqc", "run")(wildcards)
    assert actual == expected


def test_fastqc_step_part_get_output_files(ngs_data_qc):
    """Tests FastQcReportStepPart.get_output_files()"""
    # Define expected
    expected = {"fastqc_done": "work/{library_name}/report/fastqc/.done"}
    # Get actual
    actual = ngs_data_qc.get_output_files("fastqc", "run")
    assert actual == expected


def test_fastqc_step_part_get_log_file(ngs_data_qc):
    """Tests FastQcReportStepPart.get_log_file()"""
    # Define expected
    expected = "work/{library_name}/log/snakemake.fastqc.log"
    # Get actual
    actual = ngs_data_qc.get_log_file("fastqc", "run")
    assert actual == expected


def test_fastqc_step_part_get_resource_usage(ngs_data_qc):
    """Tests FastQcReportStepPart.get_resource_usage()"""
    # Define expected: default defined in workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_data_qc.get_resource("fastqc", "run", resource)
        assert actual == expected, msg_error


# Tests for NgsDataQcWorkflow ----------------------------------------------------------------------


def test_ngs_data_qc_workflow_steps(ngs_data_qc):
    """Tests simple functionality of the workflow: checks if sub steps are created."""
    # Check created sub steps
    expected = ["fastqc", "link_in", "link_out"]
    actual = list(sorted(ngs_data_qc.sub_steps.keys()))
    assert actual == expected


def test_ngs_data_qc_workflow_files(ngs_data_qc):
    """Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools."""
    # Check result file construction
    expected = [f"output/P00{i}-N1-DNA1-WGS1/report/fastqc/.done" for i in range(1, 7)]
    actual = sorted(ngs_data_qc.get_result_files())
    assert actual == expected
