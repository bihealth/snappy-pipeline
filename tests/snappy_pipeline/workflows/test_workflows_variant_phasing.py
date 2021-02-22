# -*- coding: utf-8 -*-
"""Tests for the variant_phasing workflow module code"""


import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow

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
    )


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


# Tests for PhaseByTransmissionStepPart --------------------------------------------------------------


def test_gatk_phase_by_transmission_step_part_get_input_files(variant_phasing_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_phase_by_transmission", "run")(
        wildcards
    )
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "vcf": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz",
        "vcf_md5": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "tbi": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "tbi_md5": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
    }
    assert expected == actual


def test_gatk_phase_by_transmission_step_part_get_output_files(variant_phasing_workflow):
    expected = {
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.vcf.gz.md5",
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.vcf.gz.tbi.md5",
    }
    assert (
        variant_phasing_workflow.get_output_files("gatk_phase_by_transmission", "run") == expected
    )


def test_gatk_phase_by_transmission_step_part_get_log_file(variant_phasing_workflow):
    expected = {
        "conda_info": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.conda_list.txt.md5",
        "log": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.log",
        "log_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}.log.md5",
    }
    assert variant_phasing_workflow.get_log_file("gatk_phase_by_transmission", "run") == expected


def test_gatk_phase_by_transmission_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_phase_by_transmission_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for ReadBackedPhasingOnlyStepPart --------------------------------------------------------------


def test_gatk_read_backed_phasing_only_step_part_get_input_files(variant_phasing_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_read_backed_phasing_only", "run")(
        wildcards
    )
    expected = {
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
        "tbi": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "tbi_md5": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "vcf": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz",
        "vcf_md5": "VARIANT_ANNOTATION/output/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1.vcf.gz.md5",
    }
    assert expected == actual


def test_gatk_read_backed_phasing_only_step_part_get_output_files(variant_phasing_workflow):
    expected = {
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.vcf.gz.md5",
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.vcf.gz.tbi.md5",
    }
    assert (
        variant_phasing_workflow.get_output_files("gatk_read_backed_phasing_only", "run")
        == expected
    )


def test_gatk_read_backed_phasing_only_step_part_get_log_file(variant_phasing_workflow):
    expected = {
        "conda_info": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.conda_list.txt.md5",
        "log": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.log",
        "log_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}.log.md5",
    }
    assert variant_phasing_workflow.get_log_file("gatk_read_backed_phasing_only", "run") == expected


def test_gatk_read_backed_phasing_only_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_read_backed_phasing_only_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for ReadBackedPhasingAlsoStepPart --------------------------------------------------------------


def test_gatk_read_backed_phasing_also_step_part_get_input_files(variant_phasing_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_phasing_workflow.get_input_files("gatk_read_backed_phasing_also", "run")(
        wildcards
    )
    expected = {
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
        "tbi": "work/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "tbi_md5": "work/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "vcf": "work/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1.vcf.gz",
        "vcf_md5": "work/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.P001-N1-DNA1-WGS1.vcf.gz.md5",
    }
    assert expected == actual


def test_gatk_read_backed_phasing_also_step_part_get_output_files(variant_phasing_workflow):
    expected = {
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.vcf.gz.md5",
    }
    assert (
        variant_phasing_workflow.get_output_files("gatk_read_backed_phasing_also", "run")
        == expected
    )


def test_gatk_read_backed_phasing_also_step_part_get_log_file(variant_phasing_workflow):
    expected = {
        "conda_info": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.conda_list.txt.md5",
        "log": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.log",
        "log_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}.log.md5",
    }
    assert variant_phasing_workflow.get_log_file("gatk_read_backed_phasing_also", "run") == expected


def test_gatk_read_backed_phasing_also_step_part_update_cluster_config(
    variant_phasing_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_phasing_gatk_read_backed_phasing_also_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for VariantPhasingWorkflow ----------------------------------------------------------------


def test_variant_phasing_workflow(variant_phasing_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
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
    expected = [
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.bcftools.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.freebayes.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_hc_gvcf.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.gatk_ug.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1/out/bwa.platypus.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
    ]
    expected = list(sorted(expected))
    actual = list(sorted(variant_phasing_workflow.get_result_files()))
    assert expected == actual
