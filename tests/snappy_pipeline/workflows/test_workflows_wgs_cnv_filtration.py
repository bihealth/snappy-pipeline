# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_filtration workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_filtration import WgsCnvFiltrationWorkflow

from .common import get_expected_output_vcf_files_dict
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
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          wgs_cnv_filtration:
            path_wgs_cnv_annotation: ../wgs_cnv_annotation
            tools_ngs_mapping: ['bwa']
            tools_wgs_cnv_calling: ['erds_sv2']
            # Testing 1 out 20 possible combinations:
            # {thresholds}.{inheritance}.{regions}
            filter_combinations:
              - conservative.dominant.whole_genome

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
def wgs_cnv_filtration_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return WgsCnvFiltrationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "wgs_cnv_annotation": lambda x: "WGS_CNV_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return WgsCnvFiltrationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FilterQualityStepPart ------------------------------------------------------------------


def test_filter_quality_step_part_get_input_files(wgs_cnv_filtration_workflow):
    """Tests FilterQualityStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
        }
    )
    # Define expected
    var_base_name = (
        "WGS_CNV_ANNOTATION/output/bwa.gatk3_hc.annotated.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk3_hc.annotated.P001-N1-DNA1-WGS1"
    )
    pedigree_dict = {"ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=var_base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = wgs_cnv_filtration_workflow.get_input_files("filter_quality", "run")(wildcards)
    assert actual == expected


def test_filter_quality_step_part_get_output_files(wgs_cnv_filtration_workflow):
    """Tests FilterQualityStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}/out/{mapper}.{caller}.annotated.filtered."
        r"{index_library}.{thresholds}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_filtration_workflow.get_output_files("filter_quality", "run")
    assert actual == expected


def test_filter_quality_step_part_get_log_file(wgs_cnv_filtration_workflow):
    """Tests FilterQualityStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}.{thresholds,[^\.]+}/"
        r"out/{mapper}.{caller}.annotated.filtered.{index_library}.{thresholds}.log"
    )
    actual = wgs_cnv_filtration_workflow.get_log_file("filter_quality", "run")
    assert actual == expected


def test_filter_quality_step_part_get_resource_usage(wgs_cnv_filtration_workflow):
    """Tests FilterQualityStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_filtration_workflow.get_resource("filter_quality", "run", resource)
        assert actual == expected, msg_error


# Tests for FilterInheritanceStepPart --------------------------------------------------------------


def test_filter_inheritance_step_part_get_input_files(wgs_cnv_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.annotated.filtered.P001-N1-DNA1-WGS1.conservative/out/"
        "bwa.gatk3_hc.annotated.filtered.P001-N1-DNA1-WGS1.conservative"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = wgs_cnv_filtration_workflow.get_input_files("filter_inheritance", "run")(wildcards)
    assert actual == expected


def test_filter_inheritance_step_part_get_output_files(wgs_cnv_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}/out/{mapper}.{caller}.annotated.filtered."
        r"{index_library}.{thresholds}.{inheritance}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_filtration_workflow.get_output_files("filter_inheritance", "run")
    assert actual == expected


def test_filter_inheritance_step_part_get_log_file(wgs_cnv_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}.{thresholds,[^\.]+}."
        r"{inheritance,[^\.]+}/out/{mapper}.{caller}.annotated.filtered.{index_library}."
        r"{thresholds}.{inheritance}.log"
    )
    actual = wgs_cnv_filtration_workflow.get_log_file("filter_inheritance", "run")
    assert actual == expected


def test_filter_inheritance_step_part_get_resource_usage(wgs_cnv_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_filtration_workflow.get_resource("filter_inheritance", "run", resource)
        assert actual == expected, msg_error


# Tests for FilterRegionsStepPart   ----------------------------------------------------------------


def test_filter_regions_step_part_get_input_files(wgs_cnv_filtration_workflow):
    """Tests FilterRegionsStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
            "inheritance": "dominant",
            "regions": "whole_genome",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.annotated.filtered.P001-N1-DNA1-WGS1.conservative.dominant/out/"
        "bwa.gatk3_hc.annotated.filtered.P001-N1-DNA1-WGS1.conservative.dominant"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = wgs_cnv_filtration_workflow.get_input_files("filter_regions", "run")(wildcards)
    assert actual == expected


def test_filter_regions_step_part_get_output_files(wgs_cnv_filtration_workflow):
    """Tests FilterRegionsStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}.{thresholds,[^\.]+}."
        r"{inheritance,[^\.]+}.{regions,[^\.]+}/out/{mapper}.{caller}.annotated.filtered."
        r"{index_library}.{thresholds}.{inheritance}.{regions}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_cnv_filtration_workflow.get_output_files("filter_regions", "run")
    assert actual == expected


def test_filter_regions_step_part_get_log_file(wgs_cnv_filtration_workflow):
    """Tests FilterRegionsStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}.{thresholds,[^\.]+}."
        r"{inheritance,[^\.]+}.{regions,[^\.]+}/out/{mapper}.{caller}.annotated.filtered."
        r"{index_library}.{thresholds}.{inheritance}.{regions}.log"
    )
    actual = wgs_cnv_filtration_workflow.get_log_file("filter_regions", "run")
    assert actual == expected


def test_filter_regions_step_part_get_resource_usage(wgs_cnv_filtration_workflow):
    """Tests FilterRegionsStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_filtration_workflow.get_resource("filter_regions", "run", resource)
        assert actual == expected, msg_error


# Tests for VariantFiltrationWorkflow  -------------------------------------------------------------


def test_wgs_cnv_filtration_workflow(wgs_cnv_filtration_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "filter_inheritance",
        "filter_quality",
        "filter_regions",
        "link_out",
        "write_pedigree",
    ]
    actual = list(sorted(wgs_cnv_filtration_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/bwa.erds_sv2.annotated.filtered.P00{i}-N1-DNA1-WGS1.conservative."
        "dominant.whole_genome/out/"
        "bwa.erds_sv2.annotated.filtered.P00{i}-N1-DNA1-WGS1.conservative."
        "dominant.whole_genome.{ext}"
    )
    expected = [
        tpl.format(i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
        )
    ]
    expected = sorted(expected)
    actual = sorted(wgs_cnv_filtration_workflow.get_result_files())
    assert actual == expected
