# -*- coding: utf-8 -*-
"""Tests for the variant_annotation workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow

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
          variant_annotation:
            path_jannovar_ser: /path/to/jannovar.ser

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
def variant_annotation_workflow(
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
    return VariantAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for JannovarAnnotateVcfStepPart -----------------------------------------------------------


def test_jannovar_annotate_vcf_step_part_get_input_files(variant_annotation_workflow):
    """Tests JannovarAnnotateVcfStepPart.get_input_files()"""
    # Define expected
    base_name_out = (
        "VAR_CALLING/output/"
        "{mapper}.{var_caller}.{index_ngs_library}/out/{mapper}.{var_caller}.{index_ngs_library}"
    )
    expected_vcf = base_name_out + ".vcf.gz"
    expected_tbi = base_name_out + ".vcf.gz.tbi"
    expected_keys = {"ped", "vcf", "vcf_tbi"}
    # Get actual
    result = variant_annotation_workflow.get_input_files("jannovar", "annotate_vcf")
    # Assert if all keys present
    assert set(result.keys()) == expected_keys
    # Assert vcf and tbi
    assert result["vcf"] == expected_vcf
    assert result["vcf_tbi"] == expected_tbi


def test_jannovar_annotate_vcf_step_part_get_output_files(variant_annotation_workflow):
    """Tests JannovarAnnotateVcfStepPart.get_output_files()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_annotation_workflow.get_output_files("jannovar", "annotate_vcf")
    assert actual == expected


def test_jannovar_annotate_vcf_step_part_get_log_file(variant_annotation_workflow):
    """Tests JannovarAnnotateVcfStepPart.get_log_file()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = variant_annotation_workflow.get_log_file("jannovar", "annotate_vcf")
    assert actual == expected


def test_jannovar_annotate_step_part_get_resource_usage(variant_annotation_workflow):
    """Tests JannovarAnnotateVcfStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-03:30:00", "memory": "14336M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_annotation_workflow.get_resource("jannovar", "annotate_vcf", resource)
        assert actual == expected, msg_error


# Tests for VariantAnnotationWorkflow -------------------------------------------------------------


def test_variant_annotation_workflow(variant_annotation_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["jannovar", "link_out", "write_pedigree"]
    actual = list(sorted(variant_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.jannovar_annotate_vcf.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.jannovar_annotate_vcf.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in ("vcf.gz", "vcf.gz.tbi", "vcf.gz.md5", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("gatk_hc",)
    ]
    tpl = (
        "output/{mapper}.{var_caller}.jannovar_annotate_vcf.P00{i}-N1-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.jannovar_annotate_vcf.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
        )
        for mapper in ("bwa",)
        for var_caller in ("gatk_hc",)
    ]
    actual = variant_annotation_workflow.get_result_files()
    assert actual == expected
