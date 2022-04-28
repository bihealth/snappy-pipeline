# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_signatures workflow module code"""

import textwrap

import pytest
from ruamel import yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_signatures import SomaticVariantSignaturesWorkflow

from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    return yaml.round_trip_load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fasta

          somatic_variant_signatures:
            path_somatic_variant_calling: ../SOMATIC_VARIANT_CALLING

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
def somatic_variant_signatures_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticVariantSignaturesWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return SomaticVariantSignaturesWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for TabulateVariantsStepPart ---------------------------------------------------------------


def test_tabulate_vcf_step_part_get_input_files(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_input_files()"""
    base_name = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_name + ".vcf.gz",
        "tbi": base_name + ".vcf.gz.tbi",
    }
    actual = somatic_variant_signatures_workflow.get_input_files("tabulate_vcf", "run")
    assert actual == expected


def test_tabulate_vcf_step_part_get_output_files(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_output_files()"""
    expected = {
        "tsv": (
            "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/out/"
            "{mapper}.{var_caller}.tabulate_vcf.{tumor_library}.tsv"
        )
    }
    actual = somatic_variant_signatures_workflow.get_output_files("tabulate_vcf", "run")
    assert actual == expected


def test_tabulate_vcf_step_part_get_log_file(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/log/snakemake.tabulate_vcf.log"
    )
    actual = somatic_variant_signatures_workflow.get_log_file("tabulate_vcf", "run")
    assert actual == expected


def test_tabulate_vcf_step_part_get_params(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {"tumor_library": "P001-T1-DNA1-WGS1", "normal_library": "P001-N1-DNA1-WGS1"}
    actual = somatic_variant_signatures_workflow.get_params("tabulate_vcf", "run")(wildcards)
    assert actual == expected


def test_tabulate_vcf_step_part_get_resource_usage(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "14336M", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_signatures_workflow.get_resource("tabulate_vcf", "run", resource)
        assert actual == expected, msg_error


# Tests for DeconstructSigsStepPart ----------------------------------------------------------------


def test_deconstruct_sigs_step_part_get_input_files(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_input_files()"""
    expected = {
        "tsv": (
            "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/out/"
            "{mapper}.{var_caller}.tabulate_vcf.{tumor_library}.tsv"
        )
    }
    actual = somatic_variant_signatures_workflow.get_input_files("deconstruct_sigs", "run")
    assert actual == expected


def test_deconstruct_sigs_step_part_get_output_files(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/out/"
        "{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}"
    )
    expected = {
        "tsv": base_name_out + ".tsv",
        "pdf": base_name_out + ".pdf",
    }
    actual = somatic_variant_signatures_workflow.get_output_files("deconstruct_sigs", "run")
    assert actual == expected


def test_deconstruct_sigs_step_part_get_log_file(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/log/"
        "snakemake.deconstruct_sigs.log"
    )
    actual = somatic_variant_signatures_workflow.get_log_file("deconstruct_sigs", "run")
    assert actual == expected


def test_deconstruct_sigs_step_part_get_resource_usage(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "14336M", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_signatures_workflow.get_resource(
            "deconstruct_sigs", "run", resource
        )
        assert actual == expected, msg_error


# Tests for SomaticVariantSignaturesWorkflow -------------------------------------------------------


def test_somatic_variant_signatures_workflow(somatic_variant_signatures_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["deconstruct_sigs", "link_out", "tabulate_vcf"]
    actual = list(sorted(somatic_variant_signatures_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    m_base = "output/bwa.mutect.deconstruct_sigs."
    s_base = "output/bwa.scalpel.deconstruct_sigs."
    expected = [
        m_base + "P001-T1-DNA1-WGS1/out/bwa.mutect.deconstruct_sigs.P001-T1-DNA1-WGS1.tsv",
        m_base + "P002-T1-DNA1-WGS1/out/bwa.mutect.deconstruct_sigs.P002-T1-DNA1-WGS1.tsv",
        m_base + "P002-T2-DNA1-WGS1/out/bwa.mutect.deconstruct_sigs.P002-T2-DNA1-WGS1.tsv",
        s_base + "P001-T1-DNA1-WGS1/out/bwa.scalpel.deconstruct_sigs.P001-T1-DNA1-WGS1.tsv",
        s_base + "P002-T1-DNA1-WGS1/out/bwa.scalpel.deconstruct_sigs.P002-T1-DNA1-WGS1.tsv",
        s_base + "P002-T2-DNA1-WGS1/out/bwa.scalpel.deconstruct_sigs.P002-T2-DNA1-WGS1.tsv",
    ]
    expected = set(expected)
    actual = set(somatic_variant_signatures_workflow.get_result_files())
    assert actual == expected
