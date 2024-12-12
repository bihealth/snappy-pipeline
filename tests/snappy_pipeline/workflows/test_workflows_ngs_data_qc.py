# -*- coding: utf-8 -*-
"""Tests for the ngs_data_qc workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.ngs_data_qc import NgsDataQcWorkflow

from .common import get_expected_log_files_dict
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
          ngs_mapping:
            tools:
              dna: [bwa]
            bwa:
              path_index: /path/to/bwa/index.fasta.amb
          ngs_data_qc:
            tools: ['picard']
            picard:
              path_ngs_mapping: ../ngs_mapping
              path_to_baits: /path/to/baits
              path_to_targets: /path/to/targets
              programs:
              - CollectAlignmentSummaryMetrics
              - CollectOxoGMetrics
              - CollectHsMetrics
              - CollectWgsMetrics

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
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return NgsDataQcWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_data_qc", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Construct the workflow object
    return NgsDataQcWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for PicardStepPart -------------------------------------------------------------------


def test_picard_step_part_get_output_files(ngs_data_qc):
    """Tests PicardStepPart.get_output_files() - prepare"""
    # Define expected
    expected = {
        "baits": "work/static_data/picard/out/baits.interval_list",
        "targets": "work/static_data/picard/out/targets.interval_list",
    }
    # Get actual
    actual = ngs_data_qc.get_output_files("picard", "prepare")
    assert actual == expected


def test_picard_step_part_get_log_file(ngs_data_qc):
    """Tests PicardStepPart.get_log_file() - prepare"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/static_data/picard/log/prepare",
        extended=True,
    )
    # Get actual
    actual = ngs_data_qc.get_log_file("picard", "prepare")
    assert actual == expected


def test_picard_step_part_get_input_files(ngs_data_qc):
    """Tests PicardStepPart.get_input_files() - metrics"""
    # Define expected
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "baits": "work/static_data/picard/out/baits.interval_list",
        "targets": "work/static_data/picard/out/targets.interval_list",
        "bam": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    # Get actual and assert
    actual = ngs_data_qc.get_input_files("picard", "metrics")(wildcards)
    assert actual == expected


def test_picard_step_part_get_output_files_metrics(ngs_data_qc):
    """Tests PicardStepPart.get_output_files() - metrics"""
    # Define expected
    base_out = "work/{mapper}.{library_name}/report/picard/{mapper}.{library_name}."
    expected = {
        "CollectAlignmentSummaryMetrics": base_out
        + "CollectMultipleMetrics.alignment_summary_metrics.txt",
        "CollectOxoGMetrics": base_out + "CollectOxoGMetrics.txt",
        "CollectHsMetrics": base_out + "CollectHsMetrics.txt",
        "CollectWgsMetrics": base_out + "CollectWgsMetrics.txt",
        "CollectAlignmentSummaryMetrics_md5": base_out
        + "CollectMultipleMetrics.alignment_summary_metrics.txt.md5",
        "CollectOxoGMetrics_md5": base_out + "CollectOxoGMetrics.txt.md5",
        "CollectHsMetrics_md5": base_out + "CollectHsMetrics.txt.md5",
        "CollectWgsMetrics_md5": base_out + "CollectWgsMetrics.txt.md5",
    }
    # Get actual
    actual = ngs_data_qc.get_output_files("picard", "metrics")
    assert actual == expected


def test_picard_step_part_get_log_file_metrics(ngs_data_qc):
    """Tests PicardStepPart.get_log_file() - metrics"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.{library_name}/log/picard/{mapper}.{library_name}",
        extended=True,
    )
    # Get actual
    actual = dict(ngs_data_qc.get_log_file("picard", "metrics"))
    assert actual == expected


def test_picard_step_part_get_params(ngs_data_qc):
    """Tests PicardStepPart.get_params() - metrics"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = {"prefix": "bwa.P001-N1-DNA1-WGS1"}
    # Get actual
    actual = ngs_data_qc.get_params("picard", "metrics")(wildcards)
    assert actual == expected


def test_picard_step_part_get_resource_usage(ngs_data_qc):
    """Tests PicardStepPart.get_resource_usage() - metrics"""
    # Define expected: default defined in workflow.abstract
    expected_dict = {"threads": 1, "time": "24:00:00", "memory": "24G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_data_qc.get_resource("picard", "metrics", resource)()
        assert actual == expected, msg_error


# Tests for NgsDataQcWorkflow ----------------------------------------------------------------------


def test_ngs_data_qc_workflow_steps(ngs_data_qc):
    """Tests simple functionality of the workflow: checks if sub steps are created."""
    # Check created sub steps
    expected = ["fastqc", "link_in", "link_out", "picard"]
    actual = list(sorted(ngs_data_qc.sub_steps.keys()))
    assert actual == expected


def test_ngs_data_qc_workflow_files(ngs_data_qc):
    """Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools."""
    # Check result file construction
    expected = [
        f"output/bwa.P00{i}-N1-DNA1-WGS1/report/picard/bwa.P00{i}-N1-DNA1-WGS1.{metric}.txt{ext}"
        for i in range(1, 7)
        for metric in (
            "CollectHsMetrics",
            "CollectMultipleMetrics.alignment_summary_metrics",
            "CollectOxoGMetrics",
            "CollectWgsMetrics",
        )
        for ext in ("", ".md5")
    ]
    actual = sorted(ngs_data_qc.get_result_files())
    assert actual == expected
