# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow

from .common import get_expected_log_files_dict
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
          cosmic:
            path: /path/to/cosmic.vcf.gz
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

          somatic_variant_calling:
            tools:
            - mutect
            - scalpel
            - mutect2
            scalpel:
              path_target_regions: /path/to/target/regions.bed

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
def somatic_variant_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticVariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticVariantCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for MutectStepPart ------------------------------------------------------------------------


def test_mutect_step_part_get_input_files(somatic_variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_variant_calling_workflow.get_input_files("mutect", "run")(wildcards)
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_mutect_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests file structure associated with `mutect` run in the Somatic Variant Calling Workflow."""
    # Define expected
    base_name_out = "work/{mapper}.mutect.{tumor_library}/out/{mapper}.mutect.{tumor_library}"
    expected = {
        "tbi": base_name_out + ".vcf.gz.tbi",
        "tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        "vcf": base_name_out + ".vcf.gz",
        "vcf_md5": base_name_out + ".vcf.gz.md5",
        "full_tbi": base_name_out + ".full.vcf.gz.tbi",
        "full_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        "full": base_name_out + ".full.vcf.gz",
        "full_md5": base_name_out + ".full.vcf.gz.md5",
        "txt": base_name_out + ".txt",
        "txt_md5": base_name_out + ".txt.md5",
        "wig": base_name_out + ".wig",
        "wig_md5": base_name_out + ".wig.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("mutect", "run")

    assert actual == expected


def test_mutect_step_part_get_log_file(somatic_variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.mutect.{tumor_library}/log/{mapper}.mutect.{tumor_library}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("mutect", "run")

    assert actual == expected


def test_mutect_step_part_update_cluster_config(
    somatic_variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_variant_calling_mutect_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for ScalpelStepPart ----------------------------------------------------------------------


def test_scalpel_step_part_get_input_files(somatic_variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_variant_calling_workflow.get_input_files("scalpel", "run")(wildcards)
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_scalpel_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests file structure associated with `scalpel` run in the Somatic Variant Calling
    Workflow."""
    # Define expected
    base_name_out = "work/{mapper}.scalpel.{tumor_library}/out/{mapper}.scalpel.{tumor_library}"
    expected = {
        "full_tbi": base_name_out + ".full.vcf.gz.tbi",
        "full_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        "full_vcf": base_name_out + ".full.vcf.gz",
        "full_vcf_md5": base_name_out + ".full.vcf.gz.md5",
        "tar": base_name_out + ".tar.gz",
        "tbi": base_name_out + ".vcf.gz.tbi",
        "tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        "vcf": base_name_out + ".vcf.gz",
        "vcf_md5": base_name_out + ".vcf.gz.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("scalpel", "run")

    assert actual == expected


def test_scalpel_step_part_get_log_file(somatic_variant_calling_workflow):
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.scalpel.{tumor_library}/log/{mapper}.scalpel.{tumor_library}"
    )
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("scalpel", "run")
    assert actual == expected


def test_scalpel_step_part_update_cluster_config(
    somatic_variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_variant_calling_scalpel_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for Mutect2StepPart ----------------------------------------------------------------------


def test_mutect2_step_part_get_input_files(somatic_variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "run")(wildcards)
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_mutect2_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests file structure associated with `mutect2` run in the Somatic Variant Calling
    Workflow."""
    # Define expected
    base_name_out = "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}"
    expected = {
        # Created in `filter` action
        # "full_tbi": base_name_out + ".full.vcf.gz.tbi",
        # "full_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        # "full_vcf": base_name_out + ".full.vcf.gz",
        # "full_vcf_md5": base_name_out + ".full.vcf.gz.md5",
        # "tbi": base_name_out + ".vcf.gz.tbi",
        # "tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        # "vcf": base_name_out + ".vcf.gz",
        # "vcf_md5": base_name_out + ".vcf.gz.md5",
        "raw": base_name_out + ".raw.vcf.gz",
        "raw_md5": base_name_out + ".raw.vcf.gz.md5",
        "raw_tbi": base_name_out + ".raw.vcf.gz.tbi",
        "raw_tbi_md5": base_name_out + ".raw.vcf.gz.tbi.md5",
        "stats": base_name_out + ".raw.vcf.stats",
        "stats_md5": base_name_out + ".raw.vcf.stats.md5",
        "f1r2": base_name_out + ".raw.f1r2_tar.tar.gz",
        "f1r2_md5": base_name_out + ".raw.f1r2_tar.tar.gz.md5",
        # Created in `contamination` action
        # "table": base_name_out + ".contamination.tbl",
        # "table_md5": base_name_out + ".contamination.tbl.md5",
        # "segments": base_name_out + ".segments.tbl",
        # "segments_md5": base_name_out + ".segments.tbl.md5",
        # Created in `pileup` action
        # "pileup_normal": base_name_out + ".normal.pileup",
        # "pileup_normal_md5": base_name_out + ".normal.pileup.md5",
        # "pileup_tumor": base_name_out + ".tumor.pileup",
        # "pileup_tumor_md5": base_name_out + ".tumor.pileup.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "run")

    assert actual == expected


def test_mutect2_step_part_get_log_file(somatic_variant_calling_workflow):
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.mutect2.{tumor_library}/log/{mapper}.mutect2.{tumor_library}"
    )
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "run")
    assert actual == expected


def test_mutect2_step_part_update_cluster_config(
    somatic_variant_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_variant_calling_mutect2_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for SomaticVariantCallingWorkflow ---------------------------------------------------------


def test_somatic_variant_calling_workflow(somatic_variant_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["link_out", "mutect", "scalpel", "mutect2"]
    assert set(expected).issubset(list(sorted(somatic_variant_calling_workflow.sub_steps.keys())))
    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
            "full.vcf.gz",
            "full.vcf.gz.md5",
            "full.vcf.gz.tbi",
            "full.vcf.gz.tbi.md5",
        )
        for mapper in ("bwa",)
        for var_caller in ("mutect", "scalpel", "mutect2")
    ]
    # add log files
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "conda_info.txt",
            "conda_list.txt",
            "log",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
            "log.md5",
        )
        for mapper in ("bwa",)
        for var_caller in ("mutect", "scalpel", "mutect2")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_calling_workflow.get_result_files()))
    assert expected == actual
