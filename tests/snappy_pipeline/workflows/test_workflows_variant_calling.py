# -*- coding: utf-8 -*-
"""Tests for the variant_calling workflow module code"""

from copy import deepcopy
import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    return yaml.round_trip_load(
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
            - bcftools
            - freebayes
            - gatk_hc
            - gatk_hc_gvcf
            - gatk_ug
            - platypus

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
def variant_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return VariantCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for BcftoolsStepPart ----------------------------------------------------------------------


def test_bcftools_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("bcftools", "run")(wildcards)
    expected = [
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_bcftools_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.bcftools.{index_library_name}/out/{mapper}.bcftools.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("bcftools", "run")
    assert actual == expected


def test_bcftools_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.bcftools.{index_library_name}/log/{mapper}.bcftools.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("bcftools", "run")
    assert actual == expected


def test_bcftools_step_part_update_cluster_config(variant_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["variant_calling_bcftools_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for FreebayesStepPart ---------------------------------------------------------------------


def test_freebayes_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("freebayes", "run")(wildcards)
    expected = [
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_freebayes_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.freebayes.{index_library_name}/out/{mapper}.freebayes.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("freebayes", "run")
    assert actual == expected


def test_freebayes_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.freebayes.{index_library_name}/log/{mapper}.freebayes.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("freebayes", "run")
    assert actual == expected


def test_freebayes_step_part_update_cluster_config(variant_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["variant_calling_freebayes_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for GatkHaplotypeCallerStepPart -----------------------------------------------------------


def test_gatk_hc_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk_hc", "run")(wildcards)
    expected = [
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_gatk_hc_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_hc.{index_library_name}/out/{mapper}.gatk_hc.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk_hc", "run")
    assert actual == expected


def test_gatk_hc_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_hc.{index_library_name}/log/{mapper}.gatk_hc.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("gatk_hc", "run")
    assert actual == expected


def test_gatk_hc_step_part_update_cluster_config(variant_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["variant_calling_gatk_hc_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for GatkUnifiedGenotyperStepPart ----------------------------------------------------------


def test_gatk_ug_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk_ug", "run")(wildcards)
    expected = [
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_gatk_ug_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_ug.{index_library_name}/out/{mapper}.gatk_ug.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk_ug", "run")
    assert actual == expected


def test_gatk_ug_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_ug.{index_library_name}/log/{mapper}.gatk_ug.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("gatk_ug", "run")
    assert actual == expected


def test_gatk_ug_step_part_update_cluster_config(variant_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["variant_calling_gatk_ug_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for GatkHaplotypeCallerGvcfStepPart -------------------------------------------------------


def test_gatk_hc_gvcf_step_part_discover_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk_hc_gvcf", "discover")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_gatk_hc_gvcf_step_part_discover_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_hc_gvcf.discover.{library_name}/out/"
        "{mapper}.gatk_hc_gvcf.discover.{library_name}.g"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk_hc_gvcf", "discover")
    assert actual == expected


def test_gatk_hc_gvcf_step_part_discover_get_log_file(variant_calling_workflow):
    expected = "work/{mapper}.gatk_hc_gvcf.discover.{library_name}/log/snakemake.log"
    assert variant_calling_workflow.get_log_file("gatk_hc_gvcf", "discover") == expected


def test_gatk_hc_gvcf_step_part_discover_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_gatk_hc_gvcf_discover"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


def test_gatk_hc_gvcf_step_part_combine_gvcf_get_input_files(variant_calling_workflow):
    # Define expected
    base_out = (
        "work/bwa.gatk_hc_gvcf.discover.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc_gvcf.discover.P00{i}-N1-DNA1-WGS1.g.vcf.gz"
    )
    expected = [base_out.format(i=i) for i in range(1, 7)]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    actual = variant_calling_workflow.get_input_files("gatk_hc_gvcf", "combine_gvcf")(wildcards)
    print(expected)
    assert actual == expected


def test_gatk_hc_gvcf_step_part_combine_gvcf_get_args(variant_calling_workflow):
    actual = variant_calling_workflow.get_args("gatk_hc_gvcf", "combine_gvcf")
    assert len(actual) == 1
    assert set(actual["genome_regions"].keys()) == {"1", "2"}


def test_gatk_hc_gvcf_step_part_genotype_pedigree_get_input_files(variant_calling_workflow):
    # Define expected
    base_out = (
        "work/bwa.gatk_hc_gvcf.discover.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc_gvcf.discover.P00{i}-N1-DNA1-WGS1.g.vcf.gz"
    )
    expected = [base_out.format(i=i) for i in range(1, 4)]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk_hc_gvcf", "genotype_pedigree")(
        wildcards
    )
    assert actual == expected


def test_gatk_hc_gvcf_step_part_genotype_pedigree_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk_hc_gvcf.{index_library_name}/out/"
        "{mapper}.gatk_hc_gvcf.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk_hc_gvcf", "genotype_pedigree")
    assert actual == expected


def test_gatk_hc_gvcf_step_part_genotype_pedigree_get_log_file(variant_calling_workflow):
    expected = "work/{mapper}.gatk_hc_gvcf.{index_library_name}/log/snakemake.log"
    assert variant_calling_workflow.get_log_file("gatk_hc_gvcf", "genotype_pedigree") == expected


def test_gatk_hc_gvcf_step_part_genotype_pedigree_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_gatk_hc_gvcf_genotype_pedigree"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


def test_gatk_hc_gvcf_step_part_genotype_cohort_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    actual = variant_calling_workflow.get_input_files("gatk_hc_gvcf", "genotype_cohort")(wildcards)
    expected = [
        "work/bwa.gatk_hc_gvcf.combine_gvcf/out/bwa.gatk_hc_gvcf.combine_gvcf.1.g.vcf.gz",
        "work/bwa.gatk_hc_gvcf.combine_gvcf/out/bwa.gatk_hc_gvcf.combine_gvcf.2.g.vcf.gz",
    ]
    assert actual == expected


def test_gatk_hc_gvcf_step_part_whole_cohort_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.gatk_hc_gvcf.whole_cohort/out/{mapper}.gatk_hc_gvcf.whole_cohort"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk_hc_gvcf", "genotype_cohort")
    assert actual == expected


def test_gatk_hc_gvcf_step_part_genotype_cohort_get_log_file(variant_calling_workflow):
    expected = "work/{mapper}.gatk_hc_gvcf.whole_cohort/log/snakemake.log"
    assert variant_calling_workflow.get_log_file("gatk_hc_gvcf", "genotype_cohort") == expected


def test_gatk_hc_gvcf_step_part_genotype_cohort_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_gatk_hc_gvcf_genotype_cohort"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for PlatypusStepPart ----------------------------------------------------------------------


def test_platypus_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("platypus", "run")(wildcards)
    expected = [
        "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_platypus_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.platypus.{index_library_name}/out/{mapper}.platypus.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("platypus", "run")
    assert actual == expected


def test_platypus_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.platypus.{index_library_name}/log/{mapper}.platypus.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("platypus", "run")
    assert actual == expected


def test_platypus_step_part_update_cluster_config(variant_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["variant_calling_platypus_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for VarscanStepPart ------------------------------------------------------------------------


def test_varscan_step_part_call_pedigree_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("varscan", "call_pedigree")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert expected == actual


def test_varscan_step_part_call_pedigree_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.varscan.{index_library_name}/out/{mapper}.varscan.{index_library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("varscan", "call_pedigree")
    assert actual == expected


def test_varscan_step_part_call_pedigree_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.varscan.{index_library_name}/log/{mapper}.varscan.{index_library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("varscan", "call_pedigree")

    assert actual == expected


def test_varscan_step_part_call_pedigree_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_varscan_call_pedigree"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


def test_varscan_step_part_call_cohort_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("varscan", "call_cohort")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P004-N1-DNA1-WGS1/out/bwa.P004-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P004-N1-DNA1-WGS1/out/bwa.P004-N1-DNA1-WGS1.bam.bai",
    ]
    assert expected == actual


def test_varscan_step_part_call_cohort_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.varscan.whole_cohort/out/{mapper}.varscan.whole_cohort"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_output_files("varscan", "call_cohort")
    assert actual == expected


def test_varscan_step_part_call_cohort_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.varscan.whole_cohort/log/{mapper}.varscan.whole_cohort"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_calling_workflow.get_log_file("varscan", "call_cohort")

    assert actual == expected


def test_varscan_step_part_call_cohort_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_varscan_call_cohort"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for BcftoolsStatsStepPart ------------------------------------------------------------------


def test_bcftools_stats_step_part_get_input_files(variant_calling_workflow):
    # Define expected
    vcf_file = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
    )
    expected = {"vcf": vcf_file}
    # Get actual
    actual = variant_calling_workflow.get_input_files("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/report/bcftools_stats/"
        "{mapper}.{var_caller}.{index_ngs_library}.{donor_ngs_library}"
    )
    expected = {
        "txt": base_name_out + ".txt",
        "txt_md5": base_name_out + ".txt.md5",
    }
    # Get actual
    actual = variant_calling_workflow.get_output_files("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    expected = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/log/bcftools_stats/"
        "{mapper}.{var_caller}.{index_ngs_library}.{donor_ngs_library}.log"
    )
    # Get actual
    actual = variant_calling_workflow.get_log_file("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_bcftools_stats_report"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for JannovarStatisticsStepPart -------------------------------------------------------------


def test_jannovar_statistics_step_part_get_input_files(variant_calling_workflow):
    # Define expected
    vcf_file = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
    )
    expected = {"vcf": vcf_file}
    # Get actual
    actual = variant_calling_workflow.get_input_files("jannovar_statistics", "run")
    assert actual == expected


def test_jannovar_statistics_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/report/jannovar_statistics/"
        "{mapper}.{var_caller}.{index_ngs_library}"
    )
    expected = {
        "report": base_name_out + ".txt",
        "report_md5": base_name_out + ".txt.md5",
    }
    # Get actual
    actual = variant_calling_workflow.get_output_files("jannovar_statistics", "run")
    assert actual == expected


def test_jannovar_statistics_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    expected = (
        "work/{mapper}.{var_caller}.{index_ngs_library}/log/"
        "jannovar_statistics/{mapper}.{var_caller}.{index_ngs_library}.log"
    )
    # Get actual
    actual = variant_calling_workflow.get_log_file("jannovar_statistics", "run")
    assert actual == expected


def test_jannovar_statistics_step_part_update_cluster_config(
    variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_calling_jannovar_statistics_report"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for VariantCallingWorkflow ----------------------------------------------------------------


def test_variant_calling_workflow(variant_calling_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = [
        "bcftools",
        "bcftools_stats",
        "freebayes",
        "gatk_hc",
        "gatk_hc_gvcf",
        "gatk_ug",
        "jannovar_statistics",
        "link_out",
        "platypus",
        "varscan",
        "write_pedigree",
    ]
    assert list(sorted(variant_calling_workflow.sub_steps.keys())) == expected
    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools",
            "freebayes",
            "gatk_hc",
            "gatk_hc_gvcf",
            "gatk_ug",
            "platypus",
        )
    ]
    base_out = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        base_out.format(i=i, ext=ext, mapper=mapper, var_caller=var_caller)
        for i in (1, 4)  # only for indices
        for ext in (
            "log",
            "conda_info.txt",
            "conda_list.txt",
            "log.md5",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
        )
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools",
            "freebayes",
            "gatk_hc",
            "gatk_hc_gvcf",
            "gatk_ug",
            "platypus",
        )
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "bcftools_stats/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.P00{t}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (1, 2), (1, 3), (4, 4), (4, 5), (4, 6))
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools",
            "freebayes",
            "gatk_hc",
            "gatk_hc_gvcf",
            "gatk_ug",
            "platypus",
        )
        for ext in ("txt", "txt.md5")
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "jannovar_statistics/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools",
            "freebayes",
            "gatk_hc",
            "gatk_hc_gvcf",
            "gatk_ug",
            "platypus",
        )
        for ext in ("txt", "txt.md5")
    ]
    tpl = "output/{mapper}.{var_caller}.whole_cohort/out/{mapper}.{var_caller}.whole_cohort.{ext}"
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, ext=ext)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("gatk_hc_gvcf",)
    ]
    tpl = (
        "output/{mapper}.{var_caller}.whole_cohort/report/"
        "jannovar_statistics/{mapper}.{var_caller}.whole_cohort.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, ext=ext)
        for mapper in ("bwa",)
        for var_caller in ("gatk_hc_gvcf",)
        for ext in ("txt", "txt.md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(variant_calling_workflow.get_result_files()))
    assert actual == expected


def test_variant_calling_custom_pedigree_field(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_trio_plus_sheet_fake_fs,
    mocker,
):
    """Tests VariantCallingWorkflow object pre-configured with germline trio plus sheet
    and custom pedigree field"""
    # Initialise variables
    index_standard_pedigree_list = ["P001-N1-DNA1-WGS1", "P004-N1-DNA1-WGS1", "P007-N1-DNA1-WGS1"]
    index_custom_field_pedigree_list = ["P001-N1-DNA1-WGS1", "P004-N1-DNA1-WGS1"]

    # Create alternative configuration file
    local_minimal_config = deepcopy(minimal_config)
    local_minimal_config["data_sets"]["first_batch"]["file"] = "sheet_trio_plus.tsv"

    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_trio_plus_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_trio_plus_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_calling", germline_trio_plus_sheet_fake_fs, mocker
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}

    # Construct the workflow object - should work, standard pedigree join
    vcw = VariantCallingWorkflow(
        dummy_workflow,
        local_minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )
    # Expects a single shortcut sheet
    assert len(vcw.shortcut_sheets) == 1
    sheet = vcw.shortcut_sheets[0]

    # P007 will be considered as a different index as the not defining pedigree based on `familyId`
    assert len(sheet.index_ngs_library_to_pedigree) == 3
    # Index name as expected
    assert all(
        [
            index in index_standard_pedigree_list
            for index in sheet.index_ngs_library_to_pedigree.keys()
        ]
    )

    # Construct the workflow object - should work, custom pedigree join
    local_minimal_config["data_sets"]["first_batch"]["pedigree_field"] = "familyId"
    vcw = VariantCallingWorkflow(
        dummy_workflow,
        local_minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )
    # Expects a single shortcut sheet
    assert len(vcw.shortcut_sheets) == 1
    sheet = vcw.shortcut_sheets[0]

    # P007 will be included in the same pedigree as P00{4,5,6}, i.e., all labelled 'family2'
    assert len(sheet.index_ngs_library_to_pedigree) == 2
    # Index name as expected
    assert all(
        [
            index in index_custom_field_pedigree_list
            for index in sheet.index_ngs_library_to_pedigree.keys()
        ]
    )

    # Construct the workflow object - should fail, custom pedigree field not defined
    with pytest.raises(Exception):
        local_minimal_config["data_sets"]["first_batch"]["pedigree_field"] = "_field_not_defined_"
        VariantCallingWorkflow(
            dummy_workflow,
            local_minimal_config,
            dummy_cluster_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
