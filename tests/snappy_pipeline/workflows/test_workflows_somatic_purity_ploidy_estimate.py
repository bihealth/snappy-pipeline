# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import copy
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_purity_ploidy_estimate import (
    SomaticPurityPloidyEstimateWorkflow,
)

from .conftest import patch_module_fs


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

          somatic_purity_ploidy_estimate:
            tools: ['ascat']
            tool_cnv_calling: cnvetti
            path_somatic_targeted_seq_cnv_calling: ../somatic_targeted_seq_cnv_calling
            ascat:
              b_af_loci: DUMMY

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
def minimal_config_copywritter(minimal_config):
    """Returns minimum configuration file with copywritter as the CNV caller."""
    minimal_config_adjusted = copy.deepcopy(minimal_config)
    minimal_config_adjusted["step_config"]["somatic_purity_ploidy_estimate"]["tool_cnv_calling"] = (
        "copywriter"
    )
    return minimal_config_adjusted


@pytest.fixture
def somatic_purity_ploidy_estimate_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticPurityPloidyEstimateWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)

    # Construct the workflow object
    return SomaticPurityPloidyEstimateWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def somatic_purity_ploidy_estimate_workflow_w_copywriter(
    dummy_workflow,
    minimal_config_copywritter,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticPurityPloidyEstimateWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)

    # Construct the workflow object
    return SomaticPurityPloidyEstimateWorkflow(
        dummy_workflow,
        minimal_config_copywritter,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for AscatStepPart --------------------------------------------------------------------------


def test_ascat_step_part_get_input_files_baf_tumor(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_baf_tumor()"""
    wildcards = Wildcards(fromdict={"tumor_library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "baf_tumor")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_baf_normal(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_baf_normal()"""
    wildcards = Wildcards(fromdict={"normal_library_name": "P001-N1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "baf_normal")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_cnv_tumor(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_cnv_tumor()"""
    wildcards = Wildcards(fromdict={"tumor_library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "cnv_tumor")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_cnv_normal(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_cnv_normal()"""
    wildcards = Wildcards(fromdict={"normal_library_name": "P001-N1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "cnv_normal")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_cnv_tumor_wes(
    somatic_purity_ploidy_estimate_workflow_w_copywriter,
):
    """Tests AscatStepPart._get_input_files_cnv_tumor_wes()"""
    wildcards = Wildcards(fromdict={"tumor_library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bins": (
            "../somatic_cnv_calling/work/bwa.copywriter.P001-T1-DNA1-WGS1/out/"
            "bwa.copywriter.P001-T1-DNA1-WGS1_bins.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow_w_copywriter.get_input_files(
        "ascat", "cnv_tumor_wes"
    )(wildcards)
    assert actual == expected


def test_ascat_step_part_get_input_files_cnv_normal_wes(
    somatic_purity_ploidy_estimate_workflow_w_copywriter,
):
    """Tests AscatStepPart._get_input_files_cnv_normal_wes()"""
    wildcards = Wildcards(fromdict={"normal_library_name": "P001-N1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bins": (
            "../somatic_cnv_calling/work/bwa.copywriter.P001-T1-DNA1-WGS1/out/"
            "bwa.copywriter.P001-T1-DNA1-WGS1_bins.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow_w_copywriter.get_input_files(
        "ascat", "cnv_normal_wes"
    )(wildcards)
    assert actual == expected


def test_ascat_step_part_get_output_files_baf_tumor(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_baf_tumor()"""
    expected = {
        "txt": (
            "work/{mapper}.ascat_baf_tumor.{tumor_library_name}/out/"
            "{mapper}.ascat_baf_tumor.{tumor_library_name}.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "baf_tumor")
    assert actual == expected


def test_ascat_step_part_get_output_files_baf_normal(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_baf_normal()"""
    expected = {
        "txt": (
            "work/{mapper}.ascat_baf_normal.{normal_library_name}/out/"
            "{mapper}.ascat_baf_normal.{normal_library_name}.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "baf_normal")
    assert actual == expected


def test_ascat_step_part_get_output_files_cnv_tumor(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_cnv_tumor()"""
    expected = {
        "txt": (
            "work/{mapper}.ascat_cnv_tumor.{tumor_library_name}/out/"
            "{mapper}.ascat_cnv_tumor.{tumor_library_name}.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "cnv_tumor")
    assert actual == expected


def test_ascat_step_part_get_output_files_cnv_normal(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_cnv_normal()"""
    expected = {
        "txt": (
            "work/{mapper}.ascat_cnv_normal.{normal_library_name}/out/"
            "{mapper}.ascat_cnv_normal.{normal_library_name}.txt"
        )
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "cnv_normal")
    assert actual == expected


def test_ascat_step_part_get_log_file(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart.get_log_file()"""
    # Set test cases
    base_out = (
        "work/{{mapper}}.ascat{action}.{library_type}/log/"
        "{{mapper}}.ascat{action}.{library_type}.log"
    )
    actions_and_type_dict = {
        "baf_tumor": ("_baf_tumor", "{tumor_library_name}"),
        "baf_normal": ("_baf_normal", "{normal_library_name}"),
        "cnv_tumor": ("_cnv_tumor", "{tumor_library_name}"),
        "cnv_normal": ("_cnv_normal", "{normal_library_name}"),
        "run_ascat": ("", "{tumor_library_name}"),
    }
    # Evaluate all actions
    for action, input_ in actions_and_type_dict.items():
        err_msg = f"Assert error for action '{action}'."
        expected = {"log": base_out.format(action=input_[0], library_type=input_[1])}
        actual = somatic_purity_ploidy_estimate_workflow.get_log_file("ascat", action)
        assert actual == expected, err_msg


def test_ascat_step_part_get_resource_usage(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart.get_resource()"""
    # All available actions
    actions = (
        "baf_tumor",
        "baf_normal",
        "cnv_tumor",
        "cnv_normal",
        "cnv_tumor_wes",
        "cnv_normal_wes",
        "run_ascat",
    )
    # Define expected
    expected_dict = {"threads": 8, "time": "2-00:00:00", "memory": "81920M", "partition": "medium"}
    # Evaluate
    for action in actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' for action '{action}'."
            actual = somatic_purity_ploidy_estimate_workflow.get_resource(
                "ascat", action, resource
            )()
            assert actual == expected, msg_error


# Tests for SomaticPurityPloidyEstimateWorkflow ----------------------------------------------------


def test_somatic_purity_ploidy_estimate_workflow(somatic_purity_ploidy_estimate_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["ascat", "link_out"]
    actual = list(sorted(somatic_purity_ploidy_estimate_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = [
        "output/bwa.ascat.P002-T2-DNA1-WGS1/out/.done",
        "output/bwa.ascat.P002-T1-DNA1-WGS1/out/.done",
        "output/bwa.ascat.P001-T1-DNA1-WGS1/out/.done",
    ]
    expected = set(expected)
    actual = set(somatic_purity_ploidy_estimate_workflow.get_result_files())
    assert actual == expected
