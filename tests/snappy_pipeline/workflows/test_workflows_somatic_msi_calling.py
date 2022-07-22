# -*- coding: utf-8 -*-
"""Tests for the somatic_msi_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_msi_calling import SomaticMsiCallingWorkflow

from .conftest import patch_module_fs

__author__ = "Clemens Messerschmidt"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
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
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fasta
          somatic_msi_calling:
            tools: ['mantis']
            path_ngs_mapping: ../ngs_mapping  # REQUIRED
            loci_bed: /path/to/hg19/loci.bed  # REQUIRED

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def somatic_msi_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticMsiCallingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticMsiCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FeatureCountsStepPart ------------------------------------------------------------------


def test_mantis_step_part_get_input_files(somatic_msi_calling_workflow):
    """Tests FeatureCountsStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_msi_calling_workflow.get_input_files("mantis", "run")(wildcards)
    assert actual == expected


def test_mantis_step_part_get_output_files(somatic_msi_calling_workflow):
    """Tests FeatureCountsStepPart.get_output_files()"""
    # Define expected
    base_name_out = "work/mantis.{mapper}.{library_name}/out/mantis.{mapper}.{library_name}_results"
    expected = {
        "result": base_name_out + ".txt",
        "status": base_name_out + ".txt.status",
    }
    # Get actual
    actual = somatic_msi_calling_workflow.get_output_files("mantis", "run")
    assert actual == expected


def test_mantis_step_part_get_log_file(somatic_msi_calling_workflow):
    """Tests FeatureCountsStepPart.get_log_file()"""
    expected = "work/mantis.{mapper}.{library_name}/log/snakemake.mantis_run.log"
    actual = somatic_msi_calling_workflow.get_log_file("mantis", "run")
    assert actual == expected


def test_mantis_step_part_get_resource_usage(somatic_msi_calling_workflow):
    """Tests FeatureCountsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 3, "time": "02:00:00", "memory": "92160M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_msi_calling_workflow.get_resource("mantis", "run", resource)
        assert actual == expected, msg_error


# Tests for SomaticMsiCallingWorkflow --------------------------------------------------------------


def test_somatic_msi_calling_workflow(somatic_msi_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "mantis"]
    assert list(sorted(somatic_msi_calling_workflow.sub_steps.keys())) == expected
    # Check result file construction
    expected = [
        "output/mantis.bwa.P001-T1-DNA1-WGS1/out/mantis.bwa.P001-T1-DNA1-WGS1_results.txt",
        "output/mantis.bwa.P001-T1-DNA1-WGS1/out/mantis.bwa.P001-T1-DNA1-WGS1_results.txt.status",
        "output/mantis.bwa.P002-T1-DNA1-WGS1/out/mantis.bwa.P002-T1-DNA1-WGS1_results.txt",
        "output/mantis.bwa.P002-T1-DNA1-WGS1/out/mantis.bwa.P002-T1-DNA1-WGS1_results.txt.status",
        "output/mantis.bwa.P002-T2-DNA1-WGS1/out/mantis.bwa.P002-T2-DNA1-WGS1_results.txt",
        "output/mantis.bwa.P002-T2-DNA1-WGS1/out/mantis.bwa.P002-T2-DNA1-WGS1_results.txt.status",
    ]
    actual = set(somatic_msi_calling_workflow.get_result_files())
    expected = set(expected)
    assert actual == expected
