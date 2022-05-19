# -*- coding: utf-8 -*-
"""Tests for the ``roh_calling`` workflow module code"""


import textwrap

import pytest
import ruamel.yaml as yaml
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
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

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
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for BcftoolsRohStepPart --------------------------------------------------------------------


def test_roh_calling_bcftools_roh_step_part_get_input_files_run(roh_calling_workflow):
    """Tests BcftoolsRohStepPart._get_input_files_run()"""
    # Define expected
    base_name_out = (
        "VAR_CALLING/output/bwa.gatk_hc.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.P001-N1-DNA1-WGS1"
    )
    expected = {
        "tbi": base_name_out + ".vcf.gz.tbi",
        "vcf": base_name_out + ".vcf.gz",
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    actual = roh_calling_workflow.get_input_files("bcftools_roh", "run")(wildcards)
    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_input_files_make_bed(roh_calling_workflow):
    """Tests BcftoolsRohStepPart._get_input_files_make_bed()"""
    # Define expected
    txt_file = (
        "work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz"
    )
    expected = {"txt": txt_file}
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    actual = roh_calling_workflow.get_input_files("bcftools_roh", "make_bed")(wildcards)

    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_input_files_link_bed(roh_calling_workflow):
    """Tests BcftoolsRohStepPart._get_input_files_link_bed()"""
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
    """Tests BcftoolsRohStepPart._get_output_files_run()"""
    # Define actual
    base_name_out = (
        "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}"
    )
    expected = {
        "txt": base_name_out + ".txt.gz",
        "txt_md5": base_name_out + ".txt.gz.md5",
    }
    # Get actual
    actual = roh_calling_workflow.get_output_files("bcftools_roh", "run")

    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_output_files_make_bed(roh_calling_workflow):
    """Tests BcftoolsRohStepPart._get_output_files_make_bed()"""
    expected = "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/.done"
    actual = roh_calling_workflow.get_output_files("bcftools_roh", "make_bed")
    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_shell_cmd_link_bed(roh_calling_workflow):
    """Tests BcftoolsRohStepPart.get_shell_cmd()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
            "donor_ngs_library": "P002-N1-DNA1-WGS1",
        }
    )
    base_name_str = (
        "work/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz"
    )
    expected = "\n".join(
        [
            f"ln -sr {base_name_str} {base_name_str};",
            f"ln -sr {base_name_str}.md5 {base_name_str}.md5;",
            f"ln -sr {base_name_str}.tbi {base_name_str}.tbi;",
            f"ln -sr {base_name_str}.tbi.md5 {base_name_str}.tbi.md5;",
        ]
    )
    actual = roh_calling_workflow.get_shell_cmd("bcftools_roh", "link_bed", wildcards)
    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_output_files_link_bed(roh_calling_workflow):
    """Tests BcftoolsRohStepPart._get_output_files_link_bed()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.bcftools_roh.{donor_ngs_library}"
    )
    expected = {
        "bed": base_name_out + ".bed.gz",
        "bed_md5": base_name_out + ".bed.gz.md5",
        "tbi": base_name_out + ".bed.gz.tbi",
        "tbi_md5": base_name_out + ".bed.gz.tbi.md5",
    }
    # Get actual
    actual = roh_calling_workflow.get_output_files("bcftools_roh", "link_bed")
    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_log_file(roh_calling_workflow):
    """Tests BcftoolsRohStepPart.get_log_file() - run"""
    # Define expected
    expected = (
        "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/log/"
        "snakemake.bcftools_roh.run.log"
    )
    # Get actual
    actual = roh_calling_workflow.get_log_file("bcftools_roh", "run")
    assert actual == expected


def test_roh_calling_bcftools_roh_step_part_get_resource_usage(roh_calling_workflow):
    """Tests BcftoolsRohStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "00:04:00", "memory": "4000M", "partition": "medium"}
    # Evaluate
    for action in ("run", "make_bed"):
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = roh_calling_workflow.get_resource("bcftools_roh", action, resource)
            assert actual == expected, msg_error


# Tests for RohCallingWorkflow ---------------------------------------------------------------------


def test_roh_calling_workflow(roh_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["bcftools_roh", "link_out"]
    assert expected == list(sorted(roh_calling_workflow.sub_steps.keys()))

    # Check result file construction
    p0001_base_out = "output/bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1/out/"
    p0004_base_out = "output/bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1/out/"
    expected = [
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.md5",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi",
        p0001_base_out + "bwa.gatk_hc.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.tbi",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.bed.gz.tbi.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.txt.gz",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P004-N1-DNA1-WGS1.txt.gz.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.tbi",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P005-N1-DNA1-WGS1.bed.gz.tbi.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.md5",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.tbi",
        p0004_base_out + "bwa.gatk_hc.bcftools_roh.P006-N1-DNA1-WGS1.bed.gz.tbi.md5",
    ]
    expected = list(sorted(expected))
    actual = list(sorted(roh_calling_workflow.get_result_files()))
    assert actual == expected
