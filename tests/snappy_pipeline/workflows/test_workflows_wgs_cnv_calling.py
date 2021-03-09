# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_calling workflow module code"""


import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_calling import WgsCnvCallingWorkflow

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
                dna:
                    - bwa
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
    dummy_cluster_config,
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
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for ErdsSv2StepPart -----------------------------------------------------------------------


def test_erds_sv2_step_part_get_input_files(wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_cnv_calling_workflow.get_input_files("erds_sv2", "call")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tbi": "VARIANT_CALLING/output/bwa.gatk_ug.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "vcf": "VARIANT_CALLING/output/bwa.gatk_ug.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.P001-N1-DNA1-WGS1.vcf.gz",
    }
    assert actual == expected


def test_erds_sv2_step_part_get_call_output_files(wgs_cnv_calling_workflow):
    expected = {
        "tbi": "work/{mapper}.erds_sv2.call.{library_name}/out/{mapper}.erds_sv2.call.{library_name}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.erds_sv2.call.{library_name}/out/{mapper}.erds_sv2.call.{library_name}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.erds_sv2.call.{library_name}/out/{mapper}.erds_sv2.call.{library_name}.vcf.gz",
        "vcf_md5": "work/{mapper}.erds_sv2.call.{library_name}/out/{mapper}.erds_sv2.call.{library_name}.vcf.gz.md5",
    }
    assert wgs_cnv_calling_workflow.get_output_files("erds_sv2", "call") == expected


def test_erds_sv2_step_part_get_log_file(wgs_cnv_calling_workflow):
    expected = "work/{mapper}.erds_sv2.call.{library_name}/log/snakemake.log"
    assert wgs_cnv_calling_workflow.get_log_file("erds_sv2", "call") == expected


def test_erds_sv2_step_part_update_cluster_config(wgs_cnv_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["wgs_cnv_calling_erds_sv2_call"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for VariantCallingWorkflow ----------------------------------------------------------------


def test_wgs_cnv_calling_workflow(wgs_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["cnvetti", "erds", "erds_sv2", "link_out", "write_pedigree"]
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
