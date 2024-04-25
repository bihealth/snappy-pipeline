# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_signatures workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_signatures import SomaticVariantSignaturesWorkflow

from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
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
            path_somatic_variant: ../SOMATIC_VARIANT_FILTRATION
            is_filtered: True

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
        "SOMATIC_VARIANT_FILTRATION/output/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter}.{region}/out/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter}.{region}"
    )
    expected = {
        "vcf": base_name + ".vcf.gz",
        "vcf_tbi": base_name + ".vcf.gz.tbi",
    }
    actual = somatic_variant_signatures_workflow.get_input_files("tabulate_vcf", "run")
    assert actual == expected


def test_tabulate_vcf_step_part_get_output_files(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_output_files()"""
    expected = {
        "tsv": (
            "work/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tabulate_vcf.{tumor_library}.{filter}.{region}/out/"
            "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tabulate_vcf.{tumor_library}.{filter}.{region}.tsv"
        )
    }
    actual = somatic_variant_signatures_workflow.get_output_files("tabulate_vcf", "run")
    assert actual == expected


def test_tabulate_vcf_step_part_get_log_file(somatic_variant_signatures_workflow):
    """Tests TabulateVariantsStepPart.get_log_file()"""
    expected = "work/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tabulate_vcf.{tumor_library}.{filter}.{region}/log/snakemake.tabulate_vcf.log"
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
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "14336M", "partition": "medium"}
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
            "work/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tabulate_vcf.{tumor_library}.{filter}.{region}/out/"
            "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tabulate_vcf.{tumor_library}.{filter}.{region}.tsv"
        )
    }
    actual = somatic_variant_signatures_workflow.get_input_files("deconstruct_sigs", "run")
    assert actual == expected


def test_deconstruct_sigs_step_part_get_output_files(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.deconstruct_sigs.{tumor_library}.{filter}.{region}/out/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.deconstruct_sigs.{tumor_library}.{filter}.{region}"
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
        "work/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.deconstruct_sigs.{tumor_library}.{filter}.{region}/log/"
        "snakemake.deconstruct_sigs.log"
    )
    actual = somatic_variant_signatures_workflow.get_log_file("deconstruct_sigs", "run")
    assert actual == expected


def test_deconstruct_sigs_step_part_get_resource_usage(somatic_variant_signatures_workflow):
    """Tests DeconstructSigsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "14336M", "partition": "medium"}
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
    name_pattern = "{mapper}.{caller}.{annotator}.dkfz_bias_filter.eb_filter.deconstruct_sigs.P00{i}-T{t}-DNA1-WGS1.{filt}.{region}"
    tpl = "output/" + name_pattern + "/out/" + name_pattern + ".tsv"
    expected = [
        tpl.format(
            mapper=mapper,
            caller=caller,
            annotator=annotator,
            i=i,
            t=t,
            filt=filt,
            region=region,
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for mapper in ("bwa",)
        for caller in ("mutect", "scalpel")
        for annotator in ("vep", "jannovar")
        for filt in (
            "no_filter",
            "dkfz_only",
            "dkfz_and_ebfilter",
            "dkfz_and_ebfilter_and_oxog",
            "dkfz_and_oxog",
        )
        for region in ("genome_wide",)
    ]
    expected = set(expected)
    actual = set(somatic_variant_signatures_workflow.get_result_files())
    assert actual == expected
