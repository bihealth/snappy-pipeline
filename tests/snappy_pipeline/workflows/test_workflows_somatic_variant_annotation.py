# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
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
            scalpel:
              path_target_regions: /path/to/target/regions.bed
        
          somatic_variant_annotation:
            tools: ["jannovar", "vep"]
            jannovar:
              path_jannovar_ser: /path/to/jannover.ser
            vep:
              path_dir_cache: /path/to/dir/cache

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
def somatic_variant_annotation_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticVariantAnnotationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return SomaticVariantAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for JannovarAnnotateSomaticVcfStepPart -----------------------------------------------------


def test_jannovar_step_part_get_input_files(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = somatic_variant_annotation_workflow.get_input_files("jannovar", "annotate_somatic_vcf")
    assert actual == expected


def test_jannovar_step_part_get_output_files(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}/out/"
        "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_output_files(
        "jannovar", "annotate_somatic_vcf"
    )
    assert actual == expected


def test_jannovar_step_part_get_log_file(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}/log/"
        "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_log_file("jannovar", "annotate_somatic_vcf")
    assert actual == expected


def test_jannovar_step_part_get_params(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {"tumor_library": "P001-T1-DNA1-WGS1", "normal_library": "P001-N1-DNA1-WGS1"}
    actual = somatic_variant_annotation_workflow.get_params("jannovar", "annotate_somatic_vcf")(
        wildcards
    )
    assert actual == expected


def test_jannovar_step_part_get_resource_usage(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_annotation_workflow.get_resource(
            "jannovar", "annotate_somatic_vcf", resource
        )
        assert actual == expected, msg_error


# Tests for VepAnnotateSomaticVcfStepPart ----------------------------------------------------------


def test_vep_step_part_get_input_files(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = somatic_variant_annotation_workflow.get_input_files("vep", "run")
    assert actual == expected


def test_vep_step_part_get_output_files(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.vep.{tumor_library}/out/"
        "{mapper}.{var_caller}.vep.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_output_files("vep", "run")
    assert actual == expected


def test_vep_step_part_get_log_file(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.vep.{tumor_library}/log/"
        "{mapper}.{var_caller}.vep.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_log_file("vep", "run")
    assert actual == expected


def test_vep_step_part_get_params(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {"tumor_library": "P001-T1-DNA1-WGS1", "normal_library": "P001-N1-DNA1-WGS1"}
    actual = somatic_variant_annotation_workflow.get_params("vep", "run")(wildcards)
    assert actual == expected


def test_vep_step_part_get_resource_usage(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "24:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_annotation_workflow.get_resource("vep", "run", resource)
        assert actual == expected, msg_error


# Tests for SomaticVariantAnnotationWorkflow -------------------------------------------------------


def test_somatic_variant_annotation_workflow(somatic_variant_annotation_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["jannovar", "link_out", "vep"]
    actual = list(sorted(somatic_variant_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.{annotator}.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{var_caller}.{annotator}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(
            mapper=mapper, var_caller=var_caller, annotator=annotator, i=i, t=t, ext=ext, dir_="out"
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect", "scalpel")
        for annotator in ("jannovar_annotate_somatic_vcf", "vep")
    ]
    expected += [
        tpl.format(
            mapper=mapper, var_caller=var_caller, annotator=annotator, i=i, t=t, ext=ext, dir_="log"
        )
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
        for var_caller in ("mutect", "scalpel")
        for annotator in ("jannovar_annotate_somatic_vcf", "vep")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_annotation_workflow.get_result_files()))
    assert expected == actual
