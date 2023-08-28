# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.tumor_mutational_burden import (
    TumorMutationalBurdenCalculationWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_json_files_dict
from .conftest import patch_module_fs


# Test tumor mutational burden calculation with vcf file from somatic variant calling step
@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config_calling():
    """Return YAML parsing result for configuration"""
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
            - mutect2
            - scalpel
            scalpel:
              path_target_regions: /path/to/target/regions.bed

          tumor_mutational_burden:
            path_somatic_variant_calling: ../somatic_variant_calling
            tools_ngs_mapping: []
            annotation_file: 'FALSE' # REQUIRED
            tools_somatic_variant_calling: []
            target_regions: /path/to/regions.bed

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
def tumor_mutational_burden_workflow_calling(
    dummy_workflow,
    minimal_config_calling,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return TumorMutationalBurdenCalculationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return TumorMutationalBurdenCalculationWorkflow(
        dummy_workflow,
        minimal_config_calling,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for TumorMutationalBurdenCalculationStepPart -----------------------------------------------------


def test_tumor_mutational_step_part_get_input_files_calling(tumor_mutational_burden_workflow_calling):
    """Test TumorMutationalBurdenCalculationStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = tumor_mutational_burden_workflow_calling.get_input_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_output_files_calling(tumor_mutational_burden_workflow_calling):
    """Tests TumorMutationalBurdenCalculationStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.tmb.{tumor_library}/out/"
        "{mapper}.{var_caller}.tmb.{tumor_library}"
    )
    expected = get_expected_output_json_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow_calling.get_output_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_log_files_calling(tumor_mutational_burden_workflow_calling):
    """Tests TumorMutationalBurdenCalculationStepPart.get_log_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.tmb.{tumor_library}/log/"
        "{mapper}.{var_caller}.tmb.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow_calling.get_log_file("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_resource_usage_calling(
    tumor_mutational_burden_workflow_calling,
):
    """Tests TumorMutationalBurdenCalculationStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1:00:00", "memory": "4096M"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = tumor_mutational_burden_workflow_calling.get_resource("tmb_gathering", "run", resource)
        assert actual == expected, msg_error


# Tests for TumorMutationalBurdenCalculationWorkflow -------------------------------------------------------


def test_tumor_mutational_burden_workflow_calling(tumor_mutational_burden_workflow_calling):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "tmb_gathering"]
    actual = list(sorted(tumor_mutational_burden_workflow_calling.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.tmb.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{var_caller}.tmb.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext, dir_="out")
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("json", "json.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect2", "scalpel")
    ]
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext, dir_="log")
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
        for var_caller in ("mutect2", "scalpel")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(tumor_mutational_burden_workflow_calling.get_result_files()))
    assert expected == actual


# Test Tumor mutatinal burden calculation with vcf file from somatic variant annotation step
@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config_annotation():
    """Return YAML parsing result for configuration"""
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
            - mutect2
            - scalpel
            scalpel:
              path_target_regions: /path/to/target/regions.bed

          somatic_variant_annotation:
            tools: ["jannovar", "vep"]
            jannovar:
              path_jannovar_ser: /path/to/jannover.ser
            vep:
              path_dir_cache: /path/to/dir/cache

          tumor_mutational_burden:
            path_somatic_variant_calling: ../somatic_variant_calling   # REQUIRED
            annotation_file: 'TRUE' # REQUIRED
            path_somatic_variant_annotation: ../somatic_variant_annotation  # REQUIRED
            tools_ngs_mapping: []      # default to those configured for ngs_mapping
            tools_somatic_variant_calling: []  # default to those configured for somatic_variant_calling
            tools_somatic_variant_annotation: [] # default to those configured for somatic_variant_annotation
            target_regions: /path/to/regions.bed # REQUIRED

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
def tumor_mutational_burden_workflow_annotation(
    dummy_workflow,
    minimal_config_annotation,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return TumorMutationalBurdenCalculationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
        "somatic_variant_annotation": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return TumorMutationalBurdenCalculationWorkflow(
        dummy_workflow,
        minimal_config_annotation,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for TumorMutationalBurdenCalculationStepPart -----------------------------------------------------


def test_tumor_mutational_step_part_get_input_files_annotation(tumor_mutational_burden_workflow_annotation):
    """Test TumorMutationalBurdenCalculationStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_ANNOTATION/output/{mapper}.{var_caller}.{anno_tool}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_tool}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = tumor_mutational_burden_workflow_annotation.get_input_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_output_files_annotation(tumor_mutational_burden_workflow_annotation):
    """Tests TumorMutationalBurdenCalculationStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}"
    )
    expected = get_expected_output_json_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow_annotation.get_output_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_log_files_annotation(tumor_mutational_burden_workflow_annotation):
    """Tests TumorMutationalBurdenCalculationStepPart.get_log_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}/log/"
        "{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow_annotation.get_log_file("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_resource_usage_annotation(
    tumor_mutational_burden_workflow_annotation,
):
    """Tests TumorMutationalBurdenCalculationStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1:00:00", "memory": "4096M"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = tumor_mutational_burden_workflow_annotation.get_resource("tmb_gathering", "run", resource)
        assert actual == expected, msg_error


# Tests for TumorMutationalBurdenCalculationWorkflow -------------------------------------------------------


def test_tumor_mutational_burden_workflow_annotation(tumor_mutational_burden_workflow_annotation):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "tmb_gathering"]
    actual = list(sorted(tumor_mutational_burden_workflow_annotation.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.{anno_tool}.tmb.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{var_caller}.{anno_tool}.tmb.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, anno_tool=anno_tool, i=i, t=t, ext=ext, dir_="out")
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("json", "json.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect2", "scalpel")
        for anno_tool in ("jannovar", "vep")
    ]
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, anno_tool=anno_tool, i=i, t=t, ext=ext, dir_="log")
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
        for var_caller in ("mutect2", "scalpel")
        for anno_tool in ("jannovar", "vep")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(tumor_mutational_burden_workflow_annotation.get_result_files()))
    assert expected == actual
