# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.workflows.wgs_cnv_calling import WgsCnvCallingWorkflow

from .common import (
    get_expected_gcnv_log_file,
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

          wgs_cnv_calling:
            variant_calling_tool: gatk_ug
            tools:
            - delly2
            - gcnv
            gcnv:
              precomputed_model_paths:
                - library: "default"
                  contig_ploidy: /path/to/ploidy-model
                  model_pattern: "/data/model_*"

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
    germline_sheet_fake_fs2_gcnv_model,
    mocker,
):
    """Return WgsCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs(
        "snappy_pipeline.workflows.gcnv.gcnv_run",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob.glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.gcnv.gcnv_run.glob.glob",
        return_value=["/data/model_01", "/data/model_02", "/data/model_03"],
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
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
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    base_name = (
        "work/bwa.delly2.call.P00{i}-N1-DNA1-WGS1/out/bwa.delly2.call.P00{i}-N1-DNA1-WGS1.bcf"
    )
    expected = [base_name.format(i=i) for i in (1, 2, 3)]
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "merge_calls")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_merge_calls(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - merge_calls"""
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.merge_calls.{index_ngs_library,[^\.]+}/out/"
        r"{mapper}.delly2.merge_calls.{index_ngs_library}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "merge_calls")
    assert actual == expected


def test_delly2_step_part_get_log_file_merge_calls(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_calls"""
    base_name = (
        "work/{mapper}.delly2.merge_calls.{index_ngs_library}/log/"
        "{mapper}.delly2.merge_calls.{index_ngs_library}"
    )
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
        "bcf": (
            "work/bwa.delly2.merge_calls.P001-N1-DNA1-WGS1/out/"
            "bwa.delly2.merge_calls.P001-N1-DNA1-WGS1.bcf"
        ),
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
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    base_name = (
        "work/bwa.delly2.genotype.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.delly2.genotype.P00{i}-N1-DNA1-WGS1.bcf"
    )
    expected = [base_name.format(i=i) for i in (1, 2, 3)]
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "merge_genotypes")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_merge_genotypes()"""
    base_name_out = (
        r"work/{mapper}.delly2.merge_genotypes.{index_ngs_library}/out/"
        r"{mapper}.delly2.merge_genotypes.{index_ngs_library}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "merge_genotypes")
    assert actual == expected


def test_delly2_step_part_get_log_file_merge_genotypes(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_genotypes"""
    base_name = (
        "work/{mapper}.delly2.merge_genotypes.{index_ngs_library}/log/"
        "{mapper}.delly2.merge_genotypes.{index_ngs_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "merge_genotypes")
    assert actual == expected


# Tests for Delly2StepPart (filter) ------------------


def test_delly2_step_part_get_input_files_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_filter()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    base_name = (
        "work/bwa.delly2.merge_genotypes.P001-N1-DNA1-WGS1/out/"
        "bwa.delly2.merge_genotypes.P001-N1-DNA1-WGS1"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "filter")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_output_files_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_output_files_filter()"""
    base_name_out = (
        r"work/{mapper}.delly2.filter.{index_ngs_library}/out/"
        r"{mapper}.delly2.filter.{index_ngs_library}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "filter")
    assert actual == expected


def test_delly2_step_part_get_log_file_filter(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - filter"""
    base_name = (
        "work/{mapper}.delly2.filter.{index_ngs_library}/log/"
        "{mapper}.delly2.filter.{index_ngs_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "filter")
    assert actual == expected


# Tests for Delly2StepPart (bcf_to_vcf) ------------------


def test_delly2_step_part_bcf_to_vcf_get_input_files(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_bcf_to_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    base_name = "work/bwa.delly2.filter.P001-N1-DNA1-WGS1/out/bwa.delly2.filter.P001-N1-DNA1-WGS1"
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_input_files("delly2", "bcf_to_vcf")(wildcards)
    assert actual == expected


def test_delly2_step_part_bcf_to_vcf_get_output_files(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - bcf_to_vcf"""
    base_name_out = "work/{mapper}.delly2.{library_name}/out/{mapper}.delly2.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_calling_workflow.get_output_files("delly2", "bcf_to_vcf")
    assert actual == expected


def test_delly2_step_part_bcf_to_vcf_get_log_file(wgs_cnv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - bcf_to_vcf"""
    base_name = (
        "work/{mapper}.delly2.bcf_to_vcf.{library_name}/log/"
        "{mapper}.delly2.bcf_to_vcf.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_calling_workflow.get_log_file("delly2", "bcf_to_vcf")
    assert actual == expected


# Global RunGcnvWgsStepPart Tests ------------------------------------------------------------------


def test_gcnv_call_assertion(wgs_cnv_calling_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(UnsupportedActionException):
        wgs_cnv_calling_workflow.get_input_files("gcnv", "_undefined_action_")


def test_gcnv_step_part_get_resource_usage(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_resource()"""
    # Define tested actions
    high_resource_actions = (
        "call_cnvs",
        "post_germline_calls",
    )
    all_actions = wgs_cnv_calling_workflow.substep_getattr("gcnv", "actions")
    default_actions = [action for action in all_actions if action not in high_resource_actions]
    # Define expected
    high_res_expected_dict = {
        "threads": 16,
        "time": "2-00:00:00",
        "memory": "46080M",
        "partition": "medium",
    }
    default_expected_dict = {
        "threads": 1,
        "time": "04:00:00",
        "memory": "7680M",
        "partition": "medium",
    }
    # Evaluate - high resource actions
    for action in high_resource_actions:
        for resource, expected in high_res_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = wgs_cnv_calling_workflow.get_resource("gcnv", action, resource)
            assert actual == expected, msg_error

    # Evaluate - all other actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = wgs_cnv_calling_workflow.get_resource("gcnv", action, resource)
            assert actual == expected, msg_error


def test_gcnv_get_params(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_params for all actions"""
    all_actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "joint_germline_cnv_segmentation",
    )
    actions_w_params = ("call_cnvs", "contig_ploidy", "post_germline_calls")
    for action in all_actions:
        if action in actions_w_params:
            wgs_cnv_calling_workflow.get_params("gcnv", action)
        else:
            with pytest.raises(UnsupportedActionException):
                wgs_cnv_calling_workflow.get_params("gcnv", action)


def test_gcnv_validate_precomputed_model_paths_config(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.validate_model_requirements()"""
    # Initialise input
    valid_dict = {
        "library": "library",
        "contig_ploidy": "/path/to/ploidy-model",
        "model_pattern": "/path/to/model_*",
    }
    typo_dict = {
        "library_n": "library",
        "contig_ploidy": "/path/to/ploidy-model",
        "model_pattern": "/path/to/model_*",
    }
    missing_key_dict = {"model_pattern": "/path/to/model_*"}

    # Sanity check
    wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_precomputed_model_paths_config")(
        config=[valid_dict]
    )
    # Test key typo
    with pytest.raises(InvalidConfiguration):
        wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_precomputed_model_paths_config")(
            config=[valid_dict, typo_dict]
        )
    # Test key missing
    with pytest.raises(InvalidConfiguration):
        wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_precomputed_model_paths_config")(
            config=[valid_dict, missing_key_dict]
        )


def test_gcnv_validate_ploidy_model_directory(
    fake_fs, mocker, wgs_cnv_calling_workflow, ploidy_model_files
):
    """Tests RunGcnvWgsStepPart.validate_ploidy_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/data/{file_}"
    for file_ in ploidy_model_files:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.gcnv.gcnv_run", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_ploidy_model_directory")(
        "/data"
    )
    # Should return False as empty directory
    assert not wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_ploidy_model_directory")(
        "/empty"
    )
    # Should return False not a directory
    assert not wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_ploidy_model_directory")(
        "__not_a_directory__"
    )


def test_gcnv_validate_call_model_directory(
    fake_fs, mocker, wgs_cnv_calling_workflow, call_model_files
):
    """Tests RunGcnvWgsStepPart.validate_call_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/call_data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/call_data/{file_}"
    for file_ in call_model_files:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.gcnv.gcnv_run", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_call_model_directory")(
        "/call_data"
    )
    # Should return False as empty directory
    assert not wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_call_model_directory")(
        "/empty"
    )
    # Should return False not a directory
    assert not wgs_cnv_calling_workflow.substep_getattr("gcnv", "validate_call_model_directory")(
        "__not_a_directory__"
    )


def test_gcnv_get_result_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_result_files()"""
    expected = [
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml.md5",
    ]
    actual = wgs_cnv_calling_workflow.substep_getattr("gcnv", "get_result_files")()
    assert actual == expected


# Tests for RunGcnvWgsStepPart (preprocess_intervals) ----------------------------------------------


def test_gcnv_preprocess_intervals_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_input_files_preprocess_intervals()"""
    expected = {}
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "preprocess_intervals")(None)
    assert actual == expected


def test_gcnv_preprocess_intervals_step_part_get_output_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_output_files_preprocess_intervals()"""
    output_path = (
        "work/gcnv_preprocess_intervals.{library_kit}/out/"
        "gcnv_preprocess_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": output_path}
    actual = wgs_cnv_calling_workflow.get_output_files("gcnv", "preprocess_intervals")
    assert actual == expected


def test_gcnv_target_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_log_file for 'preprocess_intervals' step"""
    expected = (
        "work/gcnv_preprocess_intervals.{library_kit}/log/"
        "gcnv_preprocess_intervals.{library_kit}.log"
    )
    actual = wgs_cnv_calling_workflow.get_log_file("gcnv", "preprocess_intervals")
    assert actual == expected


# Tests for RunGcnvWgsStepPart (coverage) ----------------------------------------------------------


def test_gcnv_coverage_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_input_files_coverage()"""
    # Define expected
    interval_list_out = (
        "work/gcnv_preprocess_intervals.default/out/"
        "gcnv_preprocess_intervals.default.interval_list"
    )
    bam_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1"
    expected = {
        "interval_list": interval_list_out,
        "bam": bam_out + ".bam",
        "bai": bam_out + ".bam.bai",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "coverage")(wildcards)
    assert actual == expected


def test_gcnv_coverage_step_part_get_output_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_output_files_coverage()"""
    tsv_out = (
        "work/{mapper}.gcnv_coverage.{library_name}/out/{mapper}.gcnv_coverage.{library_name}.tsv"
    )
    expected = {"tsv": tsv_out}
    actual = wgs_cnv_calling_workflow.get_output_files("gcnv", "coverage")
    assert actual == expected


def test_gcnv_coverage_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_log_file for 'coverage' step"""
    expected = (
        "work/{mapper}.gcnv_coverage.{library_name}/log/{mapper}.gcnv_coverage.{library_name}.log"
    )
    actual = wgs_cnv_calling_workflow.get_log_file("gcnv", "coverage")
    assert actual == expected


# Tests for RunGcnvWgsStepPart (contig_ploidy) -----------------------------------------------------


def test_gcnv_contig_ploidy_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_input_files_contig_ploidy()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_kit": "default"})
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {"tsv": tsv_list_out}
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "contig_ploidy")(wildcards)
    print(expected)
    print(actual)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_output_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_output_files_contig_ploidy()"""
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
    )
    expected = {"done": done_out}
    actual = wgs_cnv_calling_workflow.get_output_files("gcnv", "contig_ploidy")
    assert actual == expected


def test_gcnv_get_params_ploidy_model(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_params_ploidy_model()"""
    # Initialise wildcard
    wildcards = Wildcards(fromdict={"library_kit": "default"})
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test large cohort - model defined in config
    expected = {"model": "/path/to/ploidy-model"}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_ploidy_model_for_library_in_config__"}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "contig_ploidy")(wildcards_fake)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_log_file for 'contig_ploidy' step"""
    expected = get_expected_gcnv_log_file(step_name="contig_ploidy")
    actual = wgs_cnv_calling_workflow.get_log_file("gcnv", "contig_ploidy")
    assert actual == expected


# Tests for RunGcnvWgsStepPart (call_cnvs) ---------------------------------------------------------


def test_gcnv_call_cnvs_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_input_files_call_cnvs()"""
    # Define expected
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = "work/bwa.gcnv_contig_ploidy.default/out/bwa.gcnv_contig_ploidy.default/.done"
    expected = {
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_kit": "default"})
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "call_cnvs")(wildcards)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_output_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_output_files_call_cnvs()"""
    done_out = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    actual = wgs_cnv_calling_workflow.get_output_files("gcnv", "call_cnvs")
    assert actual == expected


def test_gcnv_get_params_model(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_params_model()"""
    # Initialise wildcard
    wildcards_01 = Wildcards(fromdict={"library_kit": "default", "shard": "01"})
    wildcards_02 = Wildcards(fromdict={"library_kit": "default", "shard": "02"})
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test large cohort - model defined in config - shard 01
    expected = {"model": "/data/model_01"}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "call_cnvs")(wildcards_01)
    assert actual == expected
    # Test large cohort - model defined in config - shard 02
    expected = {"model": "/data/model_02"}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "call_cnvs")(wildcards_02)
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_model_for_library_in_config__"}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "call_cnvs")(wildcards_fake)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_log_file for 'call_cnvs' step"""
    expected = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}.log"
    )
    actual = wgs_cnv_calling_workflow.get_log_file("gcnv", "call_cnvs")
    assert actual == expected


# Tests for RunGcnvWgsStepPart (post_germline_calls) -----------------------------------------------


def test_gcnv_post_germline_calls_step_part_get_input_files(
    wgs_cnv_calling_workflow,
):
    """Tests RunGcnvWgsStepPart._get_input_files_post_germline_calls()"""
    # Define expected
    call_pattern = "work/bwa.gcnv_call_cnvs.default.0{i}/out/bwa.gcnv_call_cnvs.default.0{i}/.done"
    call_list_out = [call_pattern.format(i=i) for i in range(1, 4)]  # model 01 - 03
    ploidy_out = "work/bwa.gcnv_contig_ploidy.default/out/bwa.gcnv_contig_ploidy.default/.done"
    expected = {
        "calls": call_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "post_germline_calls")(wildcards)
    assert actual == expected


def test_gcnv_get_params_post_germline_calls(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_params_post_germline_calls()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {"model": ["/data/model_01", "/data/model_02", "/data/model_03"]}
    actual = wgs_cnv_calling_workflow.get_params("gcnv", "post_germline_calls")(wildcards)
    assert actual == expected


# Tests for RunGcnvWgsStepPart (joint_germline_cnv_segmentation) -----------------------------------


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_input_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_input_files_joint_germline_cnv_segmentation()"""
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.vcf.gz"
    )
    expected = {
        "interval_list": (
            "work/gcnv_preprocess_intervals.default/out/"
            "gcnv_preprocess_intervals.default.interval_list"
        ),
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "vcf": [pattern_out.format(i=i) for i in range(1, 4)],  # P001 - P003
    }
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("gcnv", "joint_germline_cnv_segmentation")(
        wildcards
    )
    assert actual == expected


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_output_files(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart._get_output_files_joint_germline_cnv_segmentation()"""
    pattern_out = "work/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    actual = wgs_cnv_calling_workflow.get_output_files("gcnv", "joint_germline_cnv_segmentation")
    assert actual == expected


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_log_file(wgs_cnv_calling_workflow):
    """Tests RunGcnvWgsStepPart.get_log_file for 'joint_germline_cnv_segmentation' step"""
    expected = get_expected_gcnv_log_file(
        step_name="joint_germline_cnv_segmentation", extended=True
    )
    actual = wgs_cnv_calling_workflow.get_log_file("gcnv", "joint_germline_cnv_segmentation")
    assert actual == expected


# Tests for WgsCnvCallingWorkflow ------------------------------------------------------------------


def test_wgs_cnv_calling_workflow(wgs_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvetti", "delly2", "gcnv", "link_out", "write_pedigree"]
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
        for cnv_caller in ("delly2")
    ]

    # Define expected for gcnv
    pattern_out = "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/out/bwa.{tool}.P00{i}-N1-DNA1-WGS1{ext}"
    expected = [
        pattern_out.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv", "delly2")
        for ext in (
            ".vcf.gz",
            ".vcf.gz.md5",
            ".vcf.gz.tbi",
            ".vcf.gz.tbi.md5",
        )
    ]
    pattern_log = (
        "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/log/"
        "bwa.{tool}.P00{i}-N1-DNA1-WGS1.joint_germline_segmentation{ext}"
    )
    expected += [
        pattern_log.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv",)
        for ext in (
            ".conda_info.txt",
            ".conda_info.txt.md5",
            ".conda_list.txt",
            ".conda_list.txt.md5",
            ".wrapper.py",
            ".wrapper.py.md5",
            ".environment.yaml",
            ".environment.yaml.md5",
        )
    ]

    actual = wgs_cnv_calling_workflow.get_result_files()
    assert sorted(expected) == sorted(actual)
