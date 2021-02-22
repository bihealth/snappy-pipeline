# -*- coding: utf-8 -*-
"""Tests for the ``roh_calling`` workflow module code"""

from unittest.mock import Mock

import pytest
import ruamel.yaml as yaml
import sys
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.roh_calling import RohCallingWorkflow

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
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa
            star:
              path_index: /path/to/star/index
          variant_calling:
            tools:
            - gatk_hc
          roh_calling:
            path_variant_callilng: ../variant_calling

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
def roh_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return RohCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.roh_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VAR_CALLING/" + x,
    }
    # Construct the workflow object
    return RohCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for BcftoolsRohStepPart --------------------------------------------------------------------


def test_roh_calling_bcftools_roh_step_part_get_input_files_run(roh_calling_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    actual = roh_calling_workflow.get_input_files("bcftools_roh", "run")(wildcards)
    expected = {
        "tbi": "VAR_CALLING/output/bwa.gatk_hc.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "vcf": "VAR_CALLING/output/bwa.gatk_hc.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.P001-N1-DNA1-WGS1.vcf.gz",
    }
    assert expected == actual


def test_roh_calling_bcftools_roh_step_part_get_input_files_make_bed(roh_calling_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    actual = roh_calling_workflow.get_input_files("bcftools_roh", "make_bed")(wildcards)
    expected = {
        "txt": "work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz"
    }
    assert expected == actual


def test_roh_calling_bcftools_roh_step_part_get_input_files_link_bed(roh_calling_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    actual = roh_calling_workflow.get_input_files("bcftools_roh", "link_bed")(wildcards)
    expected = "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/.done"
    assert expected == actual


def test_roh_calling_bcftools_roh_step_part_get_output_files_run(roh_calling_workflow):
    expected = {
        "txt": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}.txt.gz",
        "txt_md5": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}.txt.gz.md5",
    }
    assert expected == roh_calling_workflow.get_output_files("bcftools_roh", "run")


def test_roh_calling_bcftools_roh_step_part_get_output_files_make_bed(roh_calling_workflow):
    expected = "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/.done"
    assert expected == roh_calling_workflow.get_output_files("bcftools_roh", "make_bed")


def test_roh_calling_bcftools_roh_step_part_get_shell_cmd_link_bed(roh_calling_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
            "donor_ngs_library": "P002-N1-DNA1-WGS1",
        }
    )
    expected = "\n".join(
        [
            "ln -sr work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz;",
            "ln -sr work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.md5 work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.md5;",
            "ln -sr work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi;",
            "ln -sr work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi.md5 work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi.md5;",
        ]
    )
    assert expected == roh_calling_workflow.get_shell_cmd("bcftools_roh", "link_bed", wildcards)


def test_roh_calling_bcftools_roh_step_part_get_output_files_link_bed(roh_calling_workflow):
    expected = {
        "bed": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{donor_ngs_library}.bed.gz",
        "bed_md5": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{donor_ngs_library}.bed.gz.md5",
        "tbi": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{donor_ngs_library}.bed.gz.tbi",
        "tbi_md5": "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/{mapper}.{var_caller}.bcftools_roh.{donor_ngs_library}.bed.gz.tbi.md5",
    }
    assert expected == roh_calling_workflow.get_output_files("bcftools_roh", "link_bed")


def test_roh_calling_bcftools_roh_step_part_get_log_file(roh_calling_workflow):
    expected = "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/log/snakemake.bcftools_roh.run.log"
    assert expected == roh_calling_workflow.get_log_file("bcftools_roh", "run") == expected


def test_roh_calling_bcftools_roh_step_part_run_update_cluster_config(
    roh_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["roh_calling_bcftools_roh_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual
    actual = set(dummy_cluster_config["roh_calling_bcftools_roh_make_bed"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual
    assert "roh_calling_bcftools_roh_make_link_bed" not in dummy_cluster_config


# Tests for RohCallingWorkflow ---------------------------------------------------------------------


def test_roh_calling_workflow(roh_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["bcftools_roh", "link_out"]
    assert expected == list(sorted(roh_calling_workflow.sub_steps.keys()))
    # Check result file construction
    expected = [
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.tbi.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.txt.gz",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.txt.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.tbi.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.md5",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.tbi",
        "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.tbi.md5",
    ]
    expected = list(sorted(expected))
    actual = list(sorted(roh_calling_workflow.get_result_files()))
    assert expected == actual
