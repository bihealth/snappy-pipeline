# -*- coding: utf-8 -*-
"""Tests for the somatic_hla_loh_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_hla_loh_calling import SomaticHlaLohCallingWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    return yaml.round_trip_load(
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
          somatic_hla_loh_calling:
            path_ngs_mapping: ../ngs_mapping
            path_hla_typing: ../hla_typing
            path_somatic_purity_ploidy: ../somatic_purity_ploidy_estimate

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
def somatic_hla_loh_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticHlaLohCallingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "hla_typing": lambda x: "HLA_TYPING/" + x,
    }
    # Construct the workflow object
    return SomaticHlaLohCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for LohhlaStepPart ------------------------------------------------------------------


def test_lohhla_step_part_get_input_files(somatic_hla_loh_calling_workflow):
    """Tests LohhlaStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "hla": "HLA_TYPING/output/optitype.P001-N1-DNA1-WGS1/out/optitype.P001-N1-DNA1-WGS1.txt",
    }
    actual = somatic_hla_loh_calling_workflow.get_input_files("lohhla", "run")(wildcards)
    assert actual == expected


def test_lohhla_step_part_get_output_files(somatic_hla_loh_calling_workflow):
    """Tests LohhlaStepPart.get_output_files()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.{hla_caller}.lohhla.{tumor_library}/out/"
        "{mapper}.{hla_caller}.lohhla.{tumor_library}"
    )
    expected = {"done": [base_name_out + ".done"]}
    # Get actual
    actual = somatic_hla_loh_calling_workflow.get_output_files("lohhla", "run")
    assert actual == expected


def test_lohhla_step_part_get_log_file(somatic_hla_loh_calling_workflow):
    """Tests LohhlaStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.{hla_caller}.lohhla.{tumor_library}/log/"
        "{mapper}.{hla_caller}.lohhla.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_name_log)
    actual = somatic_hla_loh_calling_workflow.get_log_file("lohhla", "run")
    assert actual == expected


def test_lohhla_step_part_get_resource_usage(somatic_hla_loh_calling_workflow):
    """Tests LohhlaStepPart.get_resource()"""
    # Define expected: default defined in workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_hla_loh_calling_workflow.get_resource("lohhla", "run", resource)
        assert actual == expected, msg_error


# Tests for SomaticHlaLohCallingWorkflow  ----------------------------------------------------------


def test_somatic_hla_loh_calling_workflow(somatic_hla_loh_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "lohhla"]
    actual = list(sorted(somatic_hla_loh_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    base_out = (
        "output/bwa.optitype.lohhla.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "bwa.optitype.lohhla.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )

    expected = [
        base_out.format(i=i, t=t, dir_="out", ext="done") for i, t in ((1, 1), (2, 1), (2, 2))
    ]
    expected += [
        base_out.format(i=i, t=t, dir_="log", ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "log",
            "conda_info.txt",
            "conda_list.txt",
            "log.md5",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
        )
    ]
    expected = set(expected)
    actual = set(somatic_hla_loh_calling_workflow.get_result_files())
    assert actual == expected
