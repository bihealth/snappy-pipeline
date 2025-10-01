# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow

from .common import get_expected_log_files_dict
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
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_calling:
            tools:
              - mutect2
            mutect2:
              contamination:
                enabled: true
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


class AttrDict(dict):
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value


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
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticVariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    scattergather = AttrDict()
    scattergather.mutect2 = lambda x: [
        x.format(scatteritem="{i}-of-24".format(i=i)) for i in range(1, 25)
    ]
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "gather": scattergather,
        "scatter": scattergather,
    }
    # Construct the workflow object
    return SomaticVariantCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Mutect2StepPart ------------------------------------------------------------------------


def test_mutect2_step_part_get_input_files_scatter(
    mutect2_wildcards, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_scatter()"""
    # Define expected
    expected = {"fai": "/path/to/ref.fa.fai"}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "scatter")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_run(mutect2_wildcards, somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_input_files_run()"""
    # Define expected
    expected = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "reference": "/path/to/ref.fa",
        "region": "work/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1/mutect2par/scatter/1-of-1.region.bed",
    }
    mutect2_wildcards_with_scatteritem = Wildcards(
        fromdict=dict(mutect2_wildcards) | {"scatteritem": "1-of-1"}
    )
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "run")(
        mutect2_wildcards_with_scatteritem
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_gather(
    mutect2_wildcards, mutect2_input_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_gather()"""
    # Define expected
    expected = {
        "vcf": [
            mutect2_input_base_name + "/mutect2par/run/{i}-of-24.raw.vcf.gz".format(i=i)
            for i in range(1, 25)
        ],
        "f1r2": [
            mutect2_input_base_name + "/mutect2par/run/{i}-of-24.raw.f1r2.tar.gz".format(i=i)
            for i in range(1, 25)
        ],
        "stats": [
            mutect2_input_base_name + "/mutect2par/run/{i}-of-24.raw.vcf.stats".format(i=i)
            for i in range(1, 25)
        ],
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "gather")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_input_files_filter(
    mutect2_wildcards, mutect2_input_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart._get_input_files_filter()"""
    # Define expected
    expected = {
        "raw": mutect2_input_base_name + ".raw.vcf.gz",
        "stats": mutect2_input_base_name + ".raw.vcf.stats",
        "orientation": mutect2_input_base_name + ".raw.read_orientation_model.tar.gz",
        "table": mutect2_input_base_name + ".contamination.tbl",
        "segments": mutect2_input_base_name + ".segments.tbl",
        "reference": "/path/to/ref.fa",
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
        "reference": "/path/to/ref.fa",
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
        "reference": "/path/to/ref.fa",
        "common_variants": "/path/to/common_variants.vcf",
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
        "reference": "/path/to/ref.fa",
        "common_variants": "/path/to/common_variants.vcf",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_input_files("mutect2", "pileup_tumor")(
        mutect2_wildcards
    )
    assert actual == expected


def test_mutect2_step_part_get_output_files_scatter(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - scatter"""
    # Define expected
    expected = {
        "regions": [
            mutect2_output_base_name + "/mutect2par/scatter/{i}-of-24.region.bed".format(i=i)
            for i in range(1, 25)
        ],
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "scatter")
    assert actual == expected


def test_mutect2_step_part_get_output_files_gather(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - gather"""
    # Define expected
    expected = {
        "vcf": mutect2_output_base_name + ".raw.vcf.gz",
        "vcf_md5": mutect2_output_base_name + ".raw.vcf.gz.md5",
        "vcf_tbi": mutect2_output_base_name + ".raw.vcf.gz.tbi",
        "vcf_tbi_md5": mutect2_output_base_name + ".raw.vcf.gz.tbi.md5",
        "stats": mutect2_output_base_name + ".raw.vcf.stats",
        "stats_md5": mutect2_output_base_name + ".raw.vcf.stats.md5",
        "orientation": mutect2_output_base_name + ".raw.read_orientation_model.tar.gz",
        "orientation_md5": mutect2_output_base_name + ".raw.read_orientation_model.tar.gz.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "gather")
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


def test_mutect2_step_part_get_output_files_run(
    mutect2_output_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_output_files() - run"""
    # Define expected
    expected = {
        "vcf": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.vcf.gz",
        "vcf_md5": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.vcf.gz.md5",
        "vcf_tbi": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.vcf.gz.tbi",
        "vcf_tbi_md5": mutect2_output_base_name
        + "/mutect2par/run/{scatteritem}.raw.vcf.gz.tbi.md5",
        "stats": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.vcf.stats",
        "stats_md5": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.vcf.stats.md5",
        "f1r2": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.f1r2.tar.gz",
        "f1r2_md5": mutect2_output_base_name + "/mutect2par/run/{scatteritem}.raw.f1r2.tar.gz.md5",
    }
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_output_files("mutect2", "run")
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


def test_mutect2_step_part_get_log_file_gather(
    mutect2_log_base_name, somatic_variant_calling_workflow
):
    """Tests Mutect2StepPart.get_log_files() - gather"""
    # Define expected
    expected = get_expected_log_files_dict(base_out=mutect2_log_base_name)
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_log_file("mutect2", "gather")
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


def test_mutect2_step_part_get_args_pileup_normal(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_pileup_normal"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    # Define expected
    expected = {"normal_lib_name": "P001-N1-DNA1-WGS1", "extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "pileup_normal")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_args_pileup_tumor(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_pileup_tumor"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    # Define expected
    expected = {"tumor_lib_name": "P001-T1-DNA1-WGS1", "extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "pileup_tumor")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_args_contamination(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_contamination"""
    # Define expected
    expected = {"extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "contamination")({})
    assert actual == expected


def test_mutect2_step_part_get_args_scatter(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_scatter"""
    # Define expected
    expected = {"padding": 5000, "ignore_chroms": [], "extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "scatter")(mutect2_wildcards)
    assert actual == expected


def test_mutect2_step_part_get_args_run(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_run"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    # Define expected
    expected = {"normal_lib_name": "P001-N1-DNA1-WGS1", "extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "run")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_args_gather(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_gather"""
    # Define expected
    expected = {}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "gather")({})
    assert actual == expected


def test_mutect2_step_part_get_args_filter(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart._get_args_filter"""
    # Define expected
    expected = {"extra_arguments": [], "java_options": ""}
    # Get actual and assert
    actual = somatic_variant_calling_workflow.get_args("mutect2", "filter")({})
    assert actual == expected


def test_mutect2_step_part_get_resource_usage_run(somatic_variant_calling_workflow):
    """Tests Mutect2StepPart.get_resource() - run"""
    # Define expected
    expected_dict = {"threads": 2, "time": "5-00:00:00", "memory": "8000M", "partition": "medium"}
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


# Tests for SomaticVariantCallingWorkflow  ---------------------------------------------------------


def test_somatic_variant_calling_workflow(somatic_variant_calling_workflow):
    """Test simple functionality of the workflow"""

    callers = {"mutect2"}
    expected = {"link_out"} | callers
    # Check created sub steps
    assert expected == set(somatic_variant_calling_workflow.sub_steps.keys())
    # Check result file construction
    mappers = ("bwa",)
    matched_tpl = (
        "output/{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    output_exts = (
        "vcf.gz",
        "vcf.gz.md5",
        "vcf.gz.tbi",
        "vcf.gz.tbi.md5",
        "full.vcf.gz",
        "full.vcf.gz.md5",
        "full.vcf.gz.tbi",
        "full.vcf.gz.tbi.md5",
    )
    expected = [
        matched_tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in output_exts
        for mapper in mappers
        for var_caller in callers
    ]
    # add log files
    matched_tpl = (
        "output/{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    meta_exts = (
        "conda_info.txt",
        "conda_list.txt",
        "log",
        "conda_info.txt.md5",
        "conda_list.txt.md5",
        "log.md5",
    )
    expected += [
        matched_tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in meta_exts
        for mapper in ("bwa",)
        for var_caller in callers
    ]

    output_exts = (
        "vcf.gz",
        "vcf.gz.md5",
        "vcf.gz.tbi",
        "vcf.gz.tbi.md5",
    )

    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_calling_workflow.get_result_files()))

    assert expected == actual
