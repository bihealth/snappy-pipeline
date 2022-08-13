# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_calling workflow module code"""


import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_calling import WgsCnvCallingWorkflow

from .common import (
    get_expected_log_files_dict,
    get_expected_output_bcf_files_dict,
    get_expected_output_vcf_files_dict,
)
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def get_expected_log_files_dict_cnvetti(base_out):
    """
    :param base_out: Base path structure for log files. For example, if the expected path for
    the log is 'work/step.path/log/step.conda_info.txt', the argument should be
    'work/step.path/log/step'.
    :type base_out: str

    :return: Returns dictionary with expected path for log files based on the provided input.
    """
    # Define expected
    expected = {
        "conda_info": base_out + ".conda_info.txt",
        "conda_list": base_out + ".conda_list.txt",
        "log": base_out + ".log",
    }
    # Return
    return expected


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
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk_ug
          wgs_cnv_calling:
            variant_calling_tool: gatk_ug
            tools:
            - erds_sv2

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
    )


@pytest.fixture
def wgs_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return WgsCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CnvettiStepPart ------------------------------------------------------------------------


def test_cnvetti_step_part_get_input_files_coverage(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    ngs_mapping_path = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    expected = {
        "bai": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam",
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "coverage")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_segment(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_segment()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    base_name = (
        "work/bwa.cnvetti_coverage.P001-N1-DNA1-WGS1/out/bwa.cnvetti_coverage.P001-N1-DNA1-WGS1"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "segment")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_merge_segments(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_merge_segments()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    base_name = (
        "work/bwa.cnvetti_segment.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.cnvetti_segment.P00{i}-N1-DNA1-WGS1.segments.{ext}"
    )
    patient_ids = (1, 2, 3, 4, 5, 6)
    expected = {
        "bcf": [base_name.format(i=i, ext="bcf") for i in patient_ids],
        "bcf_md5": [base_name.format(i=i, ext="bcf.md5") for i in patient_ids],
        "csi": [base_name.format(i=i, ext="bcf.csi") for i in patient_ids],
        "csi_md5": [base_name.format(i=i, ext="bcf.csi.md5") for i in patient_ids],
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "merge_segments")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_genotype(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_genotype()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    cov_base_name = (
        "work/bwa.cnvetti_coverage.P001-N1-DNA1-WGS1/out/bwa.cnvetti_coverage.P001-N1-DNA1-WGS1"
    )
    merge_base_name = "work/bwa.cnvetti_merge_segments/out/bwa.cnvetti_merge_segments"
    expected = {
        "sites_bcf": f"{merge_base_name}.bcf",
        "sites_bcf_md5": f"{merge_base_name}.bcf.md5",
        "sites_csi": f"{merge_base_name}.bcf.csi",
        "sites_csi_md5": f"{merge_base_name}.bcf.csi.md5",
        "coverage_bcf": f"{cov_base_name}.bcf",
        "coverage_bcf_md5": f"{cov_base_name}.bcf.md5",
        "coverage_csi": f"{cov_base_name}.bcf.csi",
        "coverage_csi_md5": f"{cov_base_name}.bcf.csi.md5",
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "genotype")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_merge_genotypes()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    base_name = (
        "work/bwa.cnvetti_genotype.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.cnvetti_genotype.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    patient_ids = (1, 2, 3, 4, 5, 6)
    expected = {
        "bcf": [base_name.format(i=i, ext="bcf") for i in patient_ids],
        "bcf_md5": [base_name.format(i=i, ext="bcf.md5") for i in patient_ids],
        "csi": [base_name.format(i=i, ext="bcf.csi") for i in patient_ids],
        "csi_md5": [base_name.format(i=i, ext="bcf.csi.md5") for i in patient_ids],
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "merge_genotypes")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_reorder_vcf(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_input_files_reorder_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    base_name = "work/bwa.cnvetti_merge_genotypes/out/bwa.cnvetti_merge_genotypes"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("cnvetti", "reorder_vcf")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_output_files_coverage(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_coverage()"""
    # Define expected
    base_name = (
        "work/{mapper}.cnvetti_coverage.{ngs_library}/out/"
        "{mapper}.cnvetti_coverage.{ngs_library}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "coverage")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_segment(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_segment()"""
    # Define expected
    seg_base_name = (
        "work/{mapper}.cnvetti_segment.{ngs_library}/out/"
        "{mapper}.cnvetti_segment.{ngs_library}.segments"
    )
    win_base_name = (
        "work/{mapper}.cnvetti_segment.{ngs_library}/out/"
        "{mapper}.cnvetti_segment.{ngs_library}.windows"
    )
    expected = {
        "segments_bcf": f"{seg_base_name}.bcf",
        "segments_bcf_md5": f"{seg_base_name}.bcf.md5",
        "segments_csi": f"{seg_base_name}.bcf.csi",
        "segments_csi_md5": f"{seg_base_name}.bcf.csi.md5",
        "windows_bcf": f"{win_base_name}.bcf",
        "windows_bcf_md5": f"{win_base_name}.bcf.md5",
        "windows_csi": f"{win_base_name}.bcf.csi",
        "windows_csi_md5": f"{win_base_name}.bcf.csi.md5",
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "segment")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_merge_segments(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_merge_segments()"""
    # Define expected
    base_name = "work/{mapper}.cnvetti_merge_segments/out/{mapper}.cnvetti_merge_segments"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "merge_segments")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_genotype(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_genotype()"""
    # Define expected
    base_name = (
        "work/{mapper}.cnvetti_genotype.{ngs_library}/out/"
        "{mapper}.cnvetti_genotype.{ngs_library}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "genotype")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_merge_genotypes()"""
    # Define expected
    base_name = "work/{mapper}.cnvetti_merge_genotypes/out/{mapper}.cnvetti_merge_genotypes"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "merge_genotypes")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_reorder_vcf(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart._get_output_files_reorder_vcf()"""
    # Define expected
    base_name = "work/{mapper}.cnvetti.{ngs_library}/out/{mapper}.cnvetti.{ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("cnvetti", "reorder_vcf")
    assert actual == expected


def test_cnvetti_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart.get_log_file()"""
    # Define test input
    merge_actions = ("merge_segments", "merge_genotypes")
    all_actions = wgs_cnv_calling_workflow.substep_getattr("cnvetti", "actions")
    default_actions = [action for action in all_actions if action not in merge_actions]

    # Define base out for dynamic expected
    merge_base_out = "work/{{mapper}}.cnvetti_{action}/log/{{mapper}}.cnvetti_{action}"
    default_base_out = (
        "work/{{mapper}}.cnvetti_{action}.{{ngs_library}}/log/"
        "{{mapper}}.cnvetti_{action}.{{ngs_library}}"
    )

    # Evaluate merge actions
    for action in merge_actions:
        b_out = merge_base_out.format(action=action)
        expected = get_expected_log_files_dict_cnvetti(base_out=b_out)
        actual = wgs_cnv_calling_workflow.get_log_file("cnvetti", action)
        assert actual == expected, f"Assert failed for action '{action}'."

    # Evaluate default actions
    for action in default_actions:
        b_out = default_base_out.format(action=action)
        expected = get_expected_log_files_dict_cnvetti(base_out=b_out)
        actual = wgs_cnv_calling_workflow.get_log_file("cnvetti", action)
        assert actual == expected, f"Assert failed for action '{action}'."


def test_cnvetti_step_part_get_ped_members(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart.get_ped_members()"""
    wildcards = Wildcards(fromdict={"ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = "P001-N1-DNA1-WGS1 P002-N1-DNA1-WGS1 P003-N1-DNA1-WGS1"
    # Get actual
    actual = wgs_cnv_calling_workflow.substep_getattr("cnvetti", "get_ped_members")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_resource_usage(wgs_cnv_calling_workflow):
    """Tests CnvettiStepPart.get_resource_usage()"""
    all_actions = wgs_cnv_calling_workflow.substep_getattr("cnvetti", "actions")
    # Define expected
    expected_dict = {"threads": 1, "time": "04:00:00", "memory": "12288M", "partition": "medium"}
    # Evaluate
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' for action '{action}'."
            actual = wgs_cnv_calling_workflow.get_resource("cnvetti", action, resource)
            assert actual == expected, msg_error


# Tests for Delly2StepPart --------------------------------------------------------------------------


def test_delly2_step_part_get_input_files_call(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_call()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    ngs_mapping_path = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    expected = {
        "bai": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam",
    }
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "call")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_call(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_call()"""
    base_name = r"work/{mapper}.delly2.call.{library_name}/out/{mapper}.delly2.call.{library_name}"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "call")
    assert actual == expected


def test_delly2_step_part_get_log_file_call(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - call"""
    base_name = "work/{mapper}.delly2.call.{library_name}/log/{mapper}.delly2.call.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "call")
    assert actual == expected


# Tests for Delly2StepPart (merge_calls) ------------------


def test_delly2_step_part_get_input_files_merge_calls(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_merge_calls()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    base_name = (
        "work/bwa.delly2.call.P00{i}-N1-DNA1-WGS1/out/bwa.delly2.call.P00{i}-N1-DNA1-WGS1.bcf"
    )
    expected = [base_name.format(i=i) for i in (1, 2, 3, 4, 5, 6)]
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "merge_calls")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_merge_calls(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - merge_calls"""
    base_name_out = r"work/{mapper,[^\.]+}.delly2.merge_calls/out/{mapper}.delly2.merge_calls"
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "merge_calls")
    assert actual == expected


def test_delly2_step_part_get_log_file_merge_calls(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_calls"""
    base_name = "work/{mapper}.delly2.merge_calls/log/{mapper}.delly2.merge_calls"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "merge_calls")
    assert actual == expected


# Tests for Delly2StepPart (genotype) ------------------


def test_delly2_step_part_get_input_files_genotype(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_genotype()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "genotype")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bcf": "work/bwa.delly2.merge_calls/out/bwa.delly2.merge_calls.bcf",
    }
    assert actual == expected


def test_delly2_step_part_get_output_files_genotype(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_genotype()"""
    base_name = (
        "work/{mapper}.delly2.genotype.{library_name}/out/{mapper}.delly2.genotype.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "genotype")
    assert actual == expected


def test_delly2_step_part_get_log_file_genotype(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - genotype"""
    base_name = (
        "work/{mapper}.delly2.genotype.{library_name}/log/{mapper}.delly2.genotype.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "genotype")
    assert actual == expected


# Tests for Delly2StepPart (merge_genotypes) ------------------


def test_delly2_step_part_get_input_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_merge_genotypes()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    base_name = (
        "work/bwa.delly2.genotype.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.delly2.genotype.P00{i}-N1-DNA1-WGS1.bcf"
    )
    expected = [base_name.format(i=i) for i in (1, 2, 3, 4, 5, 6)]
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "merge_genotypes")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_merge_genotypes()"""
    base_name_out = r"work/{mapper}.delly2.merge_genotypes/out/{mapper}.delly2.merge_genotypes"
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "merge_genotypes")
    assert actual == expected


def test_delly2_step_part_get_log_file_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_genotypes"""
    base_name = "work/{mapper}.delly2.merge_genotypes/log/{mapper}.delly2.merge_genotypes"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "merge_genotypes")
    assert actual == expected


# Tests for Delly2StepPart (filter) ------------------


def test_delly2_step_part_get_input_files_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_filter()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    base_name = "work/bwa.delly2.merge_genotypes/out/bwa.delly2.merge_genotypes"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "filter")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_filter()"""
    base_name_out = r"work/{mapper}.delly2.filter/out/{mapper}.delly2.filter"
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "filter")
    assert actual == expected


def test_delly2_step_part_get_log_file_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - filter"""
    base_name = "work/{mapper}.delly2.filter/log/{mapper}.delly2.filter"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "filter")
    assert actual == expected


# Tests for Delly2StepPart (reorder_vcf) ------------------


def test_delly2_step_part_reorder_vcf_get_input_files(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_reorder_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    base_name = "work/bwa.delly2.filter/out/bwa.delly2.filter"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "reorder_vcf")(wildcards)
    assert actual == expected


def test_delly2_step_part_reorder_vcf_get_output_files(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - reorder_vcf"""
    base_name_out = "work/{mapper}.delly2.{library_name}/out/{mapper}.delly2.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "reorder_vcf")
    assert actual == expected


def test_delly2_step_part_reorder_vcf_get_log_file(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - reorder_vcf"""
    expected = (
        "work/{mapper}.delly2.reorder_vcf.{library_name}/log/"
        "{mapper}.delly2.reorder_vcf.{library_name}.snakemake.log"
    )
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "reorder_vcf")
    assert actual == expected


# Tests for ErdsStepPart ---------------------------------------------------------------------------


def test_erds_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests ErdsStepPart.get_input_files()"""
    # Define expected
    ngs_mapping_path = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    variant_calling_path = "VARIANT_CALLING/output/bwa.gatk_ug.P001-N1-DNA1-WGS1/out/"
    expected = {
        "bai": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam",
        "tbi": variant_calling_path + "bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "vcf": variant_calling_path + "bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("erds", "run")(wildcards)
    assert actual == expected


def test_erds_step_part_get_output_files_call(wgs_cnv_calling_workflow):
    """Tests ErdsStepPart.get_output_files()"""
    # Define expected
    base_file_name = "work/{mapper}.erds.{library_name}/out/{mapper}.erds.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds", "run")
    assert actual == expected


def test_erds_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests ErdsStepPart.get_log_file()"""
    expected = "work/{mapper}.erds.{library_name}/log/snakemake.wgs_cnv_calling.log"
    actual = wgs_cnv_calling_workflow.get_log_file("erds", "run")
    assert actual == expected


def test_erds_step_part_get_resource_usage(wgs_cnv_calling_workflow):
    """Tests ErdsStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "2-00:00:00", "memory": "32768M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_calling_workflow.get_resource("erds", "run", resource)
        assert actual == expected, msg_error


# Tests for ErdsSv2StepPart -----------------------------------------------------------------------


def test_erds_sv2_step_part_get_input_files_call(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_call()"""
    # Define expected
    ngs_mapping_path = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    variant_calling_path = "VARIANT_CALLING/output/bwa.gatk_ug.P001-N1-DNA1-WGS1/out/"
    expected = {
        "bai": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": ngs_mapping_path + "bwa.P001-N1-DNA1-WGS1.bam",
        "tbi": variant_calling_path + "bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "vcf": variant_calling_path + "bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "call")(wildcards)
    assert actual == expected


def test_erds_sv2_step_part_get_input_files_merge_calls(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_merge_calls()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    expected = [
        "work/bwa.erds_sv2.call.P001-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P001-N1-DNA1-WGS1.vcf.gz",
        "work/bwa.erds_sv2.call.P002-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P002-N1-DNA1-WGS1.vcf.gz",
        "work/bwa.erds_sv2.call.P003-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P003-N1-DNA1-WGS1.vcf.gz",
        "work/bwa.erds_sv2.call.P004-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P004-N1-DNA1-WGS1.vcf.gz",
        "work/bwa.erds_sv2.call.P005-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P005-N1-DNA1-WGS1.vcf.gz",
        "work/bwa.erds_sv2.call.P006-N1-DNA1-WGS1/out/bwa.erds_sv2.call.P006-N1-DNA1-WGS1.vcf.gz",
    ]
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "merge_calls")(wildcards)
    assert actual == expected


def test_erds_sv2_step_part_get_input_files_genotype(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_genotype()"""
    lb_name = "P001-N1-DNA1-WGS1"
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": lb_name})
    # Define expected
    expected = {
        "bam": f"NGS_MAPPING/output/bwa.{lb_name}/out/bwa.{lb_name}.bam",
        "vcf_cnv": "work/bwa.erds_sv2.merge_calls/out/bwa.erds_sv2.merge_calls.vcf.gz",
        "vcf_small": (
            f"VARIANT_CALLING/output/bwa.gatk_ug.{lb_name}/out/" f"bwa.gatk_ug.{lb_name}.vcf.gz"
        ),
        "ped": f"work/write_pedigree.{lb_name}/out/{lb_name}.ped",
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "genotype")(wildcards)
    assert actual == expected


def test_erds_sv2_step_part_get_input_files_info_to_format(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_info_to_format()"""
    lb_name = "P001-N1-DNA1-WGS1"
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": lb_name})
    # Define expected
    expected = [f"work/bwa.erds_sv2.genotype.{lb_name}/out/bwa.erds_sv2.genotype.{lb_name}.vcf.gz"]
    # Get actual
    actual = sorted(
        wgs_cnv_calling_workflow.get_input_files("erds_sv2", "info_to_format")(wildcards)
    )
    assert actual == expected


def test_erds_sv2_step_part_get_input_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_merge_genotypes()"""
    work_part = "work/bwa.erds_sv2.info_to_format.P00"
    out_part = "-N1-DNA1-WGS1/out/bwa.erds_sv2.info_to_format.P00"
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    expected = [
        f"{work_part}1{out_part}1-N1-DNA1-WGS1.vcf.gz",
        f"{work_part}2{out_part}2-N1-DNA1-WGS1.vcf.gz",
        f"{work_part}3{out_part}3-N1-DNA1-WGS1.vcf.gz",
        f"{work_part}4{out_part}4-N1-DNA1-WGS1.vcf.gz",
        f"{work_part}5{out_part}5-N1-DNA1-WGS1.vcf.gz",
        f"{work_part}6{out_part}6-N1-DNA1-WGS1.vcf.gz",
    ]
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "merge_genotypes")(wildcards)
    assert actual == expected


def test_erds_sv2_step_part_get_input_files_reorder_vcf(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_input_files_reorder_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    expected = {"vcf": "work/bwa.erds_sv2.merge_genotypes/out/bwa.erds_sv2.merge_genotypes.vcf.gz"}
    # Get actual
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "reorder_vcf")(wildcards)
    assert actual == expected


def test_erds_sv2_step_part_get_output_files_call(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_output_files_call()"""
    # Define expected
    base_file_name = (
        "work/{mapper}.erds_sv2.call.{library_name}/out/{mapper}.erds_sv2.call.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds_sv2", "call")
    assert actual == expected


def test_erds_sv2_step_part_get_call_output_files_merge_calls(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_output_files_merge_calls()"""
    # Define expected
    base_name = "work/{mapper}.erds_sv2.merge_calls/out/{mapper}.erds_sv2.merge_calls"
    expected = {
        "vcf": f"{base_name}.vcf.gz",
        "tbi": f"{base_name}.vcf.gz.tbi",
        "vcf_md5": f"{base_name}.vcf.gz.md5",
        "tbi_md5": f"{base_name}.vcf.gz.tbi.md5",
    }
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds_sv2", "merge_calls")
    assert actual == expected


def test_erds_sv2_step_part_get_call_output_files_genotype(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_output_files_genotype()"""
    # Define expected
    base_file_name = (
        "work/{mapper}.erds_sv2.genotype.{library_name}/out/"
        "{mapper}.erds_sv2.genotype.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds_sv2", "genotype")
    assert actual == expected


def test_erds_sv2_step_part_get_call_output_files_info_to_format(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_output_files_info_to_format()"""
    # Define expected
    base_file_name = (
        "work/{mapper}.erds_sv2.info_to_format.{library_name}/out/"
        "{mapper}.erds_sv2.info_to_format.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds_sv2", "info_to_format")
    assert actual == expected


def test_erds_sv2_step_part_get_call_output_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart._get_output_files_merge_genotypes()"""
    # Define expected
    base_file_name = "work/{mapper}.erds_sv2.merge_genotypes/out/{mapper}.erds_sv2.merge_genotypes"
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_calling_workflow.get_output_files("erds_sv2", "merge_genotypes")
    assert actual == expected


def test_erds_sv2_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart.get_log_file()"""
    expected = "work/{mapper}.erds_sv2.call.{library_name}/log/snakemake.log"
    actual = wgs_cnv_calling_workflow.get_log_file("erds_sv2", "call")
    assert actual == expected


def test_erds_sv2_step_part_get_resource_usage(wgs_cnv_calling_workflow):
    """Tests ErdsSv2StepPart.get_resource()"""
    # Define tested actions
    cheap_actions = ("info_to_format", "merge_genotypes", "merge_calls", "reorder_vcf")
    all_actions = wgs_cnv_calling_workflow.substep_getattr("erds_sv2", "actions")
    default_actions = [action for action in all_actions if action not in cheap_actions + ("call",)]
    # Define expected
    call_expected_dict = {
        "threads": 1,
        "time": "6-08:00:00",
        "memory": "40960M",
        "partition": "medium",
    }
    cheap_expected_dict = {
        "threads": 2,
        "time": "1-00:00:00",
        "memory": "7680M",
        "partition": "medium",
    }
    default_expected_dict = {
        "threads": 4,
        "time": "6-08:00:00",
        "memory": "30720M",
        "partition": "medium",
    }
    # Evaluate - call
    for resource, expected in call_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' in action 'call'."
        actual = wgs_cnv_calling_workflow.get_resource("erds_sv2", "call", resource)
        assert actual == expected, msg_error
    # Evaluate - cheap actions
    for action in cheap_actions:
        for resource, expected in cheap_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = wgs_cnv_calling_workflow.get_resource("erds_sv2", action, resource)
            assert actual == expected, msg_error
    # Evaluate - all other actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = wgs_cnv_calling_workflow.get_resource("erds_sv2", action, resource)
            assert actual == expected, msg_error


# Tests for VariantCallingWorkflow ----------------------------------------------------------------


def test_wgs_cnv_calling_workflow(wgs_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvetti", "delly2", "erds", "erds_sv2", "link_out", "write_pedigree"]
    assert list(sorted(wgs_cnv_calling_workflow.sub_steps.keys())) == expected
    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{cnv_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, ext=ext)
        for i in (1, 4)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("erds_sv2",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_cnv_calling_workflow.get_result_files()))
    assert expected == actual
