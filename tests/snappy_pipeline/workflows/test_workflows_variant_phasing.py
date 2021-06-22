# -*- coding: utf-8 -*-
"""Tests for the variant_phasing workflow module code"""


import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    config_str = textwrap.dedent(
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
            - bcftools
            - freebayes
            - gatk_hc
            - gatk_hc_gvcf
            - gatk_ug
            - platypus
          variant_annotation:
            path_jannovar_ser: /path/to/jannovar.ser
          variant_phasing:
            path_variant_annotation: ../variant_annotation

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
    yaml_ = yaml.YAML()
    return yaml_.load(config_str)


@pytest.fixture
def variant_phasing_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantPhasingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_phasing", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VAR_CALLING/" + x,
        "variant_annotation": lambda x: "VARIANT_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return VariantPhasingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def get_expected_bam_output_file():
    """
    :return: Returns dictionary with expected bam and bai files for individuals P00{1-3}.
    """
    bam_dict = {
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
        ],
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        ],
    }
    return bam_dict


# Tests for WriteTrioPedigreeStepPart --------------------------------------------------------------


def test_write_trio_pedigree_step_part_get_output_files(variant_phasing_workflow):
    expected = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
    assert expected == variant_phasing_workflow.get_output_files("write_trio_pedigree", "run")


def test_write_trio_pedigree_step_part_run(variant_phasing_workflow, fake_fs):
    # Prepare fake file system
    fake_fs.fs.create_dir("/work/write_pedigree.P001-N1-DNA1-WGS1/out")
    # Execute trio writing
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    variant_phasing_workflow.substep_dispatch(
        "write_trio_pedigree",
        "run",
        wildcards,
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
    )
    # Check results
    assert fake_fs.os.path.exists("work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped")
    expected = (
        "\n".join(
            (
                "FAM_P001\tP001-N1-DNA1-WGS1\tP002-N1-DNA1-WGS1\tP003-N1-DNA1-WGS1\t2\t2",
                "FAM_P001\tP002-N1-DNA1-WGS1\t0\t0\t1\t1",
                "FAM_P001\tP003-N1-DNA1-WGS1\t0\t0\t2\t1",
            )
        )
        + "\n"
    )
    assert (
        expected
        == fake_fs.open("work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped").read()
    )


# Tests for PhaseByTransmissionStepPart ------------------------------------------------------------


def test_gatk_phase_by_transmission_step_part_get_input_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1"
    )
    vcf_dict = get_expected_output_vcf_files_dict(base_out=base_name_out)
    ped_dict = {"ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    expected = {**vcf_dict, **ped_dict}
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_phase_by_transmission", "run")(
        wildcards
    )
    assert actual == expected


def test_gatk_phase_by_transmission_step_part_get_output_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_output_files("gatk_phase_by_transmission", "run")
    assert actual == expected


def test_gatk_phase_by_transmission_step_part_get_log_file(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/"
        "{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_log_file("gatk_phase_by_transmission", "run")
    assert actual == expected


def test_gatk_phase_by_transmission_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_phase_by_transmission_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for ReadBackedPhasingOnlyStepPart ----------------------------------------------------------


def test_gatk_read_backed_phasing_only_step_part_get_input_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1"
    )
    vcf_dict = get_expected_output_vcf_files_dict(base_out=base_name_out)
    bam_dict = get_expected_bam_output_file()
    expected = {**bam_dict, **vcf_dict}
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_read_backed_phasing_only", "run")(
        wildcards
    )
    assert actual == expected


def test_gatk_read_backed_phasing_only_step_part_get_output_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/out/"
        "{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_output_files("gatk_read_backed_phasing_only", "run")
    assert actual == expected


def test_gatk_read_backed_phasing_only_step_part_get_log_file(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/"
        "{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_log_file("gatk_read_backed_phasing_only", "run")
    assert actual == expected


def test_gatk_read_backed_phasing_only_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_read_backed_phasing_only_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for ReadBackedPhasingAlsoStepPart ----------------------------------------------------------


def test_gatk_read_backed_phasing_also_step_part_get_input_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1"
    )
    vcf_dict = get_expected_output_vcf_files_dict(base_out=base_name_out)
    bam_dict = get_expected_bam_output_file()
    expected = {**bam_dict, **vcf_dict}
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_read_backed_phasing_also", "run")(
        wildcards
    )
    assert actual == expected


def test_gatk_read_backed_phasing_also_step_part_get_output_files(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/out/"
        "{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_output_files("gatk_read_backed_phasing_also", "run")
    assert actual == expected


def test_gatk_read_backed_phasing_also_step_part_get_log_file(variant_phasing_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/"
        "{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_phasing_workflow.get_log_file("gatk_read_backed_phasing_also", "run")
    assert actual == expected


def test_gatk_read_backed_phasing_also_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_read_backed_phasing_also_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for VariantPhasingWorkflow ----------------------------------------------------------------


def test_variant_phasing_workflow(variant_phasing_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "gatk_phase_by_transmission",
        "gatk_read_backed_phasing_also",
        "gatk_read_backed_phasing_only",
        "link_out",
        "write_trio_pedigree",
    ]
    assert expected == list(sorted(variant_phasing_workflow.sub_steps.keys()))

    # Check result file construction
    base_out = (
        "output/bwa.{tool}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.{tool}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P00{i}-N1-DNA1-WGS1"
    )
    base_out_list = [
        base_out.format(tool=tool, i=i)
        for tool in ("bcftools", "freebayes", "gatk_hc", "gatk_hc_gvcf", "gatk_ug", "platypus")
        for i in ("1", "4")
    ]
    expected = []
    for bol in base_out_list:
        tmp_dict = get_expected_output_vcf_files_dict(base_out=bol)
        expected += list(tmp_dict.values())

    expected = list(sorted(expected))
    actual = list(sorted(variant_phasing_workflow.get_result_files()))
    assert expected == actual
