# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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
            mutect2:
              common_variants: /path/to/common_variants.vcf

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
def mutect2_wildcards():
    """Returns Wildcard used in Mutect2 snakemake workflow."""
    return Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})


@pytest.fixture
def mutect2_input_base_name():
    """Returns base path used in input methods, requires only file extension addition."""
    return "work/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1"


@pytest.fixture
def mutect2_output_base_name():
    """Returns base path used in input methods, requires only file extension addition."""
    return "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}"


@pytest.fixture
def mutect2_log_base_name():
    """Returns base path used in log methods, requires only file extension addition."""
    return "work/{mapper}.mutect2.{tumor_library}/log/{mapper}.mutect2.{tumor_library}"


@pytest.fixture
def somatic_variant_calling_workflow(
    dummy_workflow,
    minimal_config,
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
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for MutectStepPart -------------------------------------------------------------------------


def test_mutect_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests MutectStepPart.get_input_files()"""
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
    """Tests MutectStepPart.get_output_files()"""
    # Define expected
    base_name_out = "work/{mapper}.mutect.{tumor_library}/out/{mapper}.mutect.{tumor_library}"
    expected = {
        "vcf_tbi": base_name_out + ".vcf.gz.tbi",
        "vcf_tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        "vcf": base_name_out + ".vcf.gz",
        "vcf_md5": base_name_out + ".vcf.gz.md5",
        "full_vcf_tbi": base_name_out + ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        "full_vcf": base_name_out + ".full.vcf.gz",
        "full_vcf_md5": base_name_out + ".full.vcf.gz.md5",
        "txt": base_name_out + ".txt",
        "txt_md5": base_name_out + ".txt.md5",
        "wig": base_name_out + ".wig",
        "wig_md5": base_name_out + ".wig.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("mutect", "run")
    assert actual == expected


def test_mutect_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests MutectStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/{mapper}.mutect.{tumor_library}/log/{mapper}.mutect.{tumor_library}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("mutect", "run")
    assert actual == expected


def test_mutect_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests MutectStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "3-00:00:00", "memory": "7577M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("mutect", "run", resource)()
        assert actual == expected, msg_error


# Tests for Mutect2StepPart ------------------------------------------------------------------------


def test_mutect2_step_part_get_input_files_run(mutect2_wildcards, somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_input_files_run()"""
    # Define expected
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "run")(mutect2_wildcards)
    assert actual == expected


def test_mutect2_step_part_get_input_files_filter(
    mutect2_wildcards, mutect2_input_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_filter()"""
    # Define expected
    expected = {
        "raw": mutect2_input_base_name + ".raw.vcf.gz",
        "stats": mutect2_input_base_name + ".raw.vcf.stats",
        "f1r2": mutect2_input_base_name + ".raw.f1r2_tar.tar.gz",
        "table": mutect2_input_base_name + ".contamination.tbl",
        "segments": mutect2_input_base_name + ".segments.tbl",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "filter")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_contamination(
    mutect2_wildcards, mutect2_input_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_contamination()"""
    # Define expected
    expected = {
        "normal": mutect2_input_base_name + ".normal.pileup",
        "tumor": mutect2_input_base_name + ".tumor.pileup",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "contamination")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_pileup_normal(
    mutect2_wildcards, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_pileup_normal()"""
    # Define expected
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "pileup_normal")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_pileup_tumor(
    mutect2_wildcards, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_pileup_tumor()"""
    # Define expected
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "pileup_tumor")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_output_files_run(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - run"""
    # Define expected
    expected = {
        "raw": mutect2_output_base_name + ".raw.vcf.gz",
        "raw_md5": mutect2_output_base_name + ".raw.vcf.gz.md5",
        "raw_tbi": mutect2_output_base_name + ".raw.vcf.gz.tbi",
        "raw_tbi_md5": mutect2_output_base_name + ".raw.vcf.gz.tbi.md5",
        "stats": mutect2_output_base_name + ".raw.vcf.stats",
        "stats_md5": mutect2_output_base_name + ".raw.vcf.stats.md5",
        "f1r2": mutect2_output_base_name + ".raw.f1r2_tar.tar.gz",
        "f1r2_md5": mutect2_output_base_name + ".raw.f1r2_tar.tar.gz.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "run")
    assert actual == expected


def test_mutect2_step_part_get_output_files_filter(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - filter"""
    # Define expected
    expected = {
        "full_vcf": mutect2_output_base_name + ".full.vcf.gz",
        "full_vcf_md5": mutect2_output_base_name + ".full.vcf.gz.md5",
        "full_vcf_tbi": mutect2_output_base_name + ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": mutect2_output_base_name + ".full.vcf.gz.tbi.md5",
        "vcf": mutect2_output_base_name + ".vcf.gz",
        "vcf_md5": mutect2_output_base_name + ".vcf.gz.md5",
        "vcf_tbi": mutect2_output_base_name + ".vcf.gz.tbi",
        "vcf_tbi_md5": mutect2_output_base_name + ".vcf.gz.tbi.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "filter")
    assert actual == expected


def test_mutect2_step_part_get_output_files_contamination(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - contamination"""
    # Define expected
    expected = {
        "table": mutect2_output_base_name + ".contamination.tbl",
        "table_md5": mutect2_output_base_name + ".contamination.tbl.md5",
        "segments": mutect2_output_base_name + ".segments.tbl",
        "segments_md5": mutect2_output_base_name + ".segments.tbl.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "contamination")
    assert actual == expected


def test_mutect2_step_part_get_output_files_pileup_normal(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - pileup_normal"""
    # Define expected
    expected = {
        "pileup": mutect2_output_base_name + ".normal.pileup",
        "pileup_md5": mutect2_output_base_name + ".normal.pileup.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "pileup_normal")
    assert actual == expected


def test_mutect2_step_part_get_output_files_pileup_tumor(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - pileup_tumor"""
    # Define expected
    expected = {
        "pileup": mutect2_output_base_name + ".tumor.pileup",
        "pileup_md5": mutect2_output_base_name + ".tumor.pileup.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "pileup_tumor")
    assert actual == expected


def test_mutect2_step_part_get_log_file_run(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - run"""
    # Define expected
    expected = get_expected_log_files_dict(base_out=mutect2_log_base_name)
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "run")
    assert actual == expected


def test_mutect2_step_part_get_log_file_filter(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - filter"""
    # Define expected
    expected = {
        "log": mutect2_log_base_name + ".filter.log",
        "log_md5": mutect2_log_base_name + ".filter.log.md5",
        "conda_info": mutect2_log_base_name + ".filter.conda_info.txt",
        "conda_info_md5": mutect2_log_base_name + ".filter.conda_info.txt.md5",
        "conda_list": mutect2_log_base_name + ".filter.conda_list.txt",
        "conda_list_md5": mutect2_log_base_name + ".filter.conda_list.txt.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "filter")
    assert actual == expected


def test_mutect2_step_part_get_log_file_contamination(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - contamination"""
    # Define expected
    expected = {
        "log": mutect2_log_base_name + ".contamination.log",
        "log_md5": mutect2_log_base_name + ".contamination.log.md5",
        "conda_info": mutect2_log_base_name + ".contamination.conda_info.txt",
        "conda_info_md5": mutect2_log_base_name + ".contamination.conda_info.txt.md5",
        "conda_list": mutect2_log_base_name + ".contamination.conda_list.txt",
        "conda_list_md5": mutect2_log_base_name + ".contamination.conda_list.txt.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "contamination")
    assert actual == expected


def test_mutect2_step_part_get_log_file_pileup_normal(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - pileup_normal"""
    # Define expected
    expected = {
        "log": mutect2_log_base_name + ".pileup_normal.log",
        "log_md5": mutect2_log_base_name + ".pileup_normal.log.md5",
        "conda_info": mutect2_log_base_name + ".pileup_normal.conda_info.txt",
        "conda_info_md5": mutect2_log_base_name + ".pileup_normal.conda_info.txt.md5",
        "conda_list": mutect2_log_base_name + ".pileup_normal.conda_list.txt",
        "conda_list_md5": mutect2_log_base_name + ".pileup_normal.conda_list.txt.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "pileup_normal")
    assert actual == expected


def test_mutect2_step_part_get_log_file_pileup_tumor(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - pileup_tumor"""
    # Define expected
    expected = {
        "log": mutect2_log_base_name + ".pileup_tumor.log",
        "log_md5": mutect2_log_base_name + ".pileup_tumor.log.md5",
        "conda_info": mutect2_log_base_name + ".pileup_tumor.conda_info.txt",
        "conda_info_md5": mutect2_log_base_name + ".pileup_tumor.conda_info.txt.md5",
        "conda_list": mutect2_log_base_name + ".pileup_tumor.conda_list.txt",
        "conda_list_md5": mutect2_log_base_name + ".pileup_tumor.conda_list.txt.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "pileup_tumor")
    assert actual == expected


def test_mutect2_step_part_get_resource_usage_run(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - run"""
    # Define expected
    expected_dict = {"threads": 2, "time": "5-00:00:00", "memory": "3584M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("mutect2", "run", resource)()
        assert actual == expected, msg_error


def test_mutect2_step_part_get_resource_usage_filter(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - filter"""
    # Define expected
    expected_dict = {"threads": 2, "time": "03:59:00", "memory": "15872M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("mutect2", "filter", resource)()
        assert actual == expected, msg_error


def test_mutect2_step_part_get_resource_usage_contamination(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - contamination"""
    # Define expected
    expected_dict = {"threads": 2, "time": "03:59:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource(
            "mutect2", "contamination", resource
        )()
        assert actual == expected, msg_error


def test_mutect2_step_part_get_resource_usage_pileup_normal(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - pileup_normal"""
    # Define expected
    expected_dict = {"threads": 2, "time": "03:59:00", "memory": "8000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource(
            "mutect2", "pileup_normal", resource
        )()
        assert actual == expected, msg_error


def test_mutect2_step_part_get_resource_usage_pileup_tumor(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - pileup_tumor"""
    # Define expected
    expected_dict = {"threads": 2, "time": "03:59:00", "memory": "8000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource(
            "mutect2", "pileup_tumor", resource
        )()
        assert actual == expected, msg_error


# Tests for ScalpelStepPart ----------------------------------------------------------------------


def test_scalpel_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests ScalpelStepPart.get_input_files()"""
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
    """Tests ScalpelStepPart.get_output_files()"""
    # Define expected
    base_name_out = "work/{mapper}.scalpel.{tumor_library}/out/{mapper}.scalpel.{tumor_library}"
    expected = {
        "full_vcf_tbi": base_name_out + ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        "full_vcf": base_name_out + ".full.vcf.gz",
        "full_vcf_md5": base_name_out + ".full.vcf.gz.md5",
        "tar": base_name_out + ".tar.gz",
        "tar_md5": base_name_out + ".tar.gz.md5",
        "vcf_tbi": base_name_out + ".vcf.gz.tbi",
        "vcf_tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        "vcf": base_name_out + ".vcf.gz",
        "vcf_md5": base_name_out + ".vcf.gz.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("scalpel", "run")
    assert actual == expected


def test_scalpel_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests ScalpelStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.scalpel.{tumor_library}/log/{mapper}.scalpel.{tumor_library}"
    )
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("scalpel", "run")
    assert actual == expected


def test_scalpel_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests ScalpelStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "81920M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("scalpel", "run", resource)()
        assert actual == expected, msg_error


# Tests for Strelka2StepPart ----------------------------------------------------------------------


def test_strelka2_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests Strelka2StepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_variant_calling_workflow.get_input_files("strelka2", "run")(wildcards)
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_strelka2_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests Strelka2StepPart.get_output_files()"""
    # Define expected
    base_name_out = "work/{mapper}.strelka2.{tumor_library}/out/{mapper}.strelka2.{tumor_library}"
    expected = {
        "vcf": base_name_out + ".vcf.gz",
        "vcf_md5": base_name_out + ".vcf.gz.md5",
        "vcf_tbi": base_name_out + ".vcf.gz.tbi",
        "vcf_tbi_md5": base_name_out + ".vcf.gz.tbi.md5",
        "full_vcf": base_name_out + ".full.vcf.gz",
        "full_vcf_md5": base_name_out + ".full.vcf.gz.md5",
        "full_vcf_tbi": base_name_out + ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": base_name_out + ".full.vcf.gz.tbi.md5",
        "stats": base_name_out + ".tsv",
        "stats_md5": base_name_out + ".tsv.md5",
        "report": base_name_out + ".xml",
        "report_md5": base_name_out + ".xml.md5",
        "bed": base_name_out + ".bed.gz",
        "bed_md5": base_name_out + ".bed.gz.md5",
        "bed_tbi": base_name_out + ".bed.gz.tbi",
        "bed_tbi_md5": base_name_out + ".bed.gz.tbi.md5",
    }
    # Get actual
    actual = somatic_variant_calling_workflow.get_output_files("strelka2", "run")
    assert actual == expected


def test_strelka2_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests Strelka2StepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.strelka2.{tumor_library}/log/{mapper}.strelka2.{tumor_library}"
    )
    # Get actual
    actual = somatic_variant_calling_workflow.get_log_file("strelka2", "run")
    assert actual == expected


def test_strelka2_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests Strelka2StepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "4G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("strelka2", "run", resource)()
        assert actual == expected, msg_error


# Tests for BcftoolsJointStepPart   ----------------------------------------------------------------


def test_bcftools_joint_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests BcftoolsJointStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "donor_name": "P001"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam.bai",
        ],
    }
    actual = somatic_variant_calling_workflow.get_input_files("bcftools_joint", "run")(wildcards)
    assert actual == expected


def test_bcftools_joint_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests BcftoolsJointStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.bcftools_joint.{donor_name}/out/{mapper}.bcftools_joint.{donor_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = somatic_variant_calling_workflow.get_output_files("bcftools_joint", "run")
    assert actual == expected


def test_bcftools_joint_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests BcftoolsJointStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.bcftools_joint.{donor_name}/log/{mapper}.bcftools_joint.{donor_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_log)
    actual = somatic_variant_calling_workflow.get_log_file("bcftools_joint", "run")
    assert actual == expected


def test_bcftools_joint_step_part_get_args(somatic_variant_calling_workflow):
    """Tests BcftoolsJointStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"donor_name": "P001"})
    expected = {
        "sample_list": ["P001-N1-DNA1-WGS1", "P001-T1-DNA1-WGS1", "P001-T1-RNA1-mRNA_seq1"],
        "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
    }
    actual = somatic_variant_calling_workflow.get_args("bcftools_joint", "run")(wildcards)
    assert actual == expected


def test_bcftools_joint_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests BcftoolsJointStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("bcftools_joint", "run", resource)()
        assert actual == expected, msg_error


# Tests for VarscanJointStepPart -------------------------------------------------------------------


def test_varscan_joint_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests VarscanJointStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "donor_name": "P001"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam.bai",
        ],
    }
    actual = somatic_variant_calling_workflow.get_input_files("varscan_joint", "run")(wildcards)
    assert actual == expected


def test_varscan_joint_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests VarscanJointStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.varscan_joint.{donor_name}/out/{mapper}.varscan_joint.{donor_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = somatic_variant_calling_workflow.get_output_files("varscan_joint", "run")
    assert actual == expected


def test_varscan_joint_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests VarscanJointStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.varscan_joint.{donor_name}/log/{mapper}.varscan_joint.{donor_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_log)
    actual = somatic_variant_calling_workflow.get_log_file("varscan_joint", "run")
    assert actual == expected


def test_varscan_joint_step_part_get_args(somatic_variant_calling_workflow):
    """Tests VarscanJointStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"donor_name": "P001"})
    expected = {
        "sample_list": ["P001-N1-DNA1-WGS1", "P001-T1-DNA1-WGS1", "P001-T1-RNA1-mRNA_seq1"],
        "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
    }
    actual = somatic_variant_calling_workflow.get_args("varscan_joint", "run")(wildcards)
    assert actual == expected


def test_varscan_joint_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests VarscanJointStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "2-00:00:00", "memory": "1024M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("varscan_joint", "run", resource)()
        assert actual == expected, msg_error


# Tests for PlatypusJointStepPart  -----------------------------------------------------------------


def test_platypus_joint_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests PlatypusJointStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "donor_name": "P001"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam.bai",
        ],
    }
    actual = somatic_variant_calling_workflow.get_input_files("platypus_joint", "run")(wildcards)
    assert actual == expected


def test_platypus_joint_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests PlatypusJointStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.platypus_joint.{donor_name}/out/{mapper}.platypus_joint.{donor_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = somatic_variant_calling_workflow.get_output_files("platypus_joint", "run")
    assert actual == expected


def test_platypus_joint_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests PlatypusJointStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.platypus_joint.{donor_name}/log/{mapper}.platypus_joint.{donor_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_log)
    actual = somatic_variant_calling_workflow.get_log_file("platypus_joint", "run")
    assert actual == expected


def test_platypus_joint_step_part_get_args(somatic_variant_calling_workflow):
    """Tests PlatypusJointStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"donor_name": "P001"})
    expected = {
        "sample_list": ["P001-N1-DNA1-WGS1", "P001-T1-DNA1-WGS1", "P001-T1-RNA1-mRNA_seq1"],
        "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
    }
    actual = somatic_variant_calling_workflow.get_args("platypus_joint", "run")(wildcards)
    assert actual == expected


def test_platypus_joint_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests PlatypusJointStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "61440M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("platypus_joint", "run", resource)()
        assert actual == expected, msg_error


# Tests for GatkHcJointStepPart  -------------------------------------------------------------------


def test_gatk_hc_joint_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests GatkHcJointStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "donor_name": "P001"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam.bai",
        ],
    }
    actual = somatic_variant_calling_workflow.get_input_files("gatk_hc_joint", "run")(wildcards)
    assert actual == expected


def test_gatk_hc_joint_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests GatkHcJointStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.gatk_hc_joint.{donor_name}/out/{mapper}.gatk_hc_joint.{donor_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = somatic_variant_calling_workflow.get_output_files("gatk_hc_joint", "run")
    assert actual == expected


def test_gatk_hc_joint_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests GatkHcJointStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.gatk_hc_joint.{donor_name}/log/{mapper}.gatk_hc_joint.{donor_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_log)
    actual = somatic_variant_calling_workflow.get_log_file("gatk_hc_joint", "run")
    assert actual == expected


def test_gatk_hc_joint_step_part_get_args(somatic_variant_calling_workflow):
    """Tests GatkHcJointStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"donor_name": "P001"})
    expected = {
        "sample_list": ["P001-N1-DNA1-WGS1", "P001-T1-DNA1-WGS1", "P001-T1-RNA1-mRNA_seq1"],
        "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
    }
    actual = somatic_variant_calling_workflow.get_args("gatk_hc_joint", "run")(wildcards)
    assert actual == expected


def test_gatk_hc_joint_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests GatkHcJointStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "3-08:00:00", "memory": "2048M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("gatk_hc_joint", "run", resource)()
        assert actual == expected, msg_error


# Tests for GatkUgJointStepPart  -------------------------------------------------------------------


def test_gatk_ug_joint_step_part_get_input_files(somatic_variant_calling_workflow):
    """Tests GatkUgJointStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "donor_name": "P001"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-RNA1-mRNA_seq1/out/bwa.P001-T1-RNA1-mRNA_seq1.bam.bai",
        ],
    }
    actual = somatic_variant_calling_workflow.get_input_files("gatk_ug_joint", "run")(wildcards)
    assert actual == expected


def test_gatk_ug_joint_step_part_get_output_files(somatic_variant_calling_workflow):
    """Tests GatkUgJointStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.gatk_ug_joint.{donor_name}/out/{mapper}.gatk_ug_joint.{donor_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = somatic_variant_calling_workflow.get_output_files("gatk_ug_joint", "run")
    assert actual == expected


def test_gatk_ug_joint_step_part_get_log_file(somatic_variant_calling_workflow):
    """Tests GatkUgJointStepPart.get_log_file()"""
    base_name_log = (
        "work/{mapper}.gatk_ug_joint.{donor_name}/log/{mapper}.gatk_ug_joint.{donor_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_log)
    actual = somatic_variant_calling_workflow.get_log_file("gatk_ug_joint", "run")
    assert actual == expected


def test_gatk_ug_joint_step_part_get_args(somatic_variant_calling_workflow):
    """Tests GatkUgJointStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"donor_name": "P001"})
    expected = {
        "sample_list": ["P001-N1-DNA1-WGS1", "P001-T1-DNA1-WGS1", "P001-T1-RNA1-mRNA_seq1"],
        "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
    }
    actual = somatic_variant_calling_workflow.get_args("gatk_ug_joint", "run")(wildcards)
    assert actual == expected


def test_gatk_ug_joint_step_part_get_resource_usage(somatic_variant_calling_workflow):
    """Tests GatkUgJointStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "3-08:00:00", "memory": "2048M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_calling_workflow.get_resource("gatk_ug_joint", "run", resource)()
        assert actual == expected, msg_error


# Tests for SomaticVariantCallingWorkflow  ---------------------------------------------------------


def test_somatic_variant_calling_workflow(somatic_variant_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["link_out", "mutect", "scalpel"]
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
    # add special cases
    expected += [
        tpl.format(mapper=mapper, var_caller="mutect", i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("txt", "txt.md5", "wig", "wig.md5")
        for mapper in ("bwa",)
    ]
    expected += [
        tpl.format(mapper=mapper, var_caller="scalpel", i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("tar.gz", "tar.gz.md5")
        for mapper in ("bwa",)
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
