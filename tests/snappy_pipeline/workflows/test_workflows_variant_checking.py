# -*- coding: utf-8 -*-
"""Tests for the variant_checking workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.variant_checking import VariantCheckingWorkflow

from .conftest import patch_module_fs


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

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk3_hc
            gatk3_hc: {}

          variant_checking:
            tools_ngs_mapping: ['bwa']  # optional, copied from ngs mapping config
            tools_variant_calling: ['gatk3_hc']  # optional, copied from variant calling config
            path_variant_calling: ../variant_calling  # REQUIRED
            tools: ['peddy']

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
def variant_checking_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return VariantAnnotationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VAR_CALLING/" + x,
    }
    # Construct the workflow object
    return VariantCheckingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for PeddyStepPart --------------------------------------------------------------------------


def test_peddy_step_part_get_input_files(variant_checking_workflow):
    """Tests PeddyStepPart.get_input_files()"""
    # Define expected
    base_name_out = (
        "VAR_CALLING/output/"
        "{mapper}.{var_caller}.{index_ngs_library}/out/{mapper}.{var_caller}.{index_ngs_library}"
    )
    expected_vcf = base_name_out + ".vcf.gz"
    expected_tbi = base_name_out + ".vcf.gz.tbi"
    expected_keys = {"ped", "vcf", "vcf_tbi"}
    # Get actual
    result = variant_checking_workflow.get_input_files("peddy", "run")
    # Assert if all keys present
    assert set(result.keys()) == expected_keys
    # Assert vcf and tbi
    assert result["vcf"] == expected_vcf
    assert result["vcf_tbi"] == expected_tbi


def test_peddy_step_part_get_output_files(variant_checking_workflow):
    """Tests PeddyStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{var_caller}.peddy.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.peddy.{index_ngs_library}"
    )
    expected = {
        "background_pca": base_name_out + ".background_pca.json",
        "het_check": base_name_out + ".het_check.csv",
        "html": base_name_out + ".html",
        "ped_check": base_name_out + ".ped_check.csv",
        "ped": base_name_out + ".peddy.ped",
        "sex_check": base_name_out + ".sex_check.csv",
    }
    actual = variant_checking_workflow.get_output_files("peddy", "run")
    assert actual == expected


def test_peddy_step_part_get_log_file(variant_checking_workflow):
    """Tests PeddyStepPart.get_log_file()"""
    expected = "work/{mapper}.{var_caller}.peddy.{index_ngs_library}/log/snakemake.filter.log"
    actual = variant_checking_workflow.get_log_file("peddy", "run")
    assert actual == expected


def test_peddy_step_part_get_resource_usage(variant_checking_workflow):
    """Tests PeddyStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "10:00:00", "memory": "15360M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_checking_workflow.get_resource("peddy", "run", resource)()
        assert actual == expected, msg_error


# Tests for VariantCheckingWorkflow ----------------------------------------------------------------


def test_variant_checking_workflow(variant_checking_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "peddy", "write_pedigree"]
    actual = list(sorted(variant_checking_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "work/{mapper}.{var_caller}.peddy.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.peddy.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            "background_pca.json",
            "het_check.csv",
            "html",
            "ped_check.csv",
            "peddy.ped",
            "sex_check.csv",
        )
        for mapper in ("bwa",)
        for var_caller in ("gatk3_hc",)
    ]
    actual = variant_checking_workflow.get_result_files()
    assert actual == expected
