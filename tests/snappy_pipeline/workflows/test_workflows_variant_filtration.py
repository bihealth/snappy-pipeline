# -*- coding: utf-8 -*-
"""Tests for the variant_filtration workflow module code"""

from pathlib import Path
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


def get_project_root():
    """Get project root

    :return: Return path to project root as a string.
    """
    return str(Path(__file__).parent.parent.parent.parent)


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
            - gatk3_hc

          variant_filtration:
            path_variant_annotation: ../variant_annotation
            tools_ngs_mapping: ['bwa']
            tools_variant_calling: ['gatk3_hc']
            # Testing 1 out 40+ possible combinations:
            # {thresholds}.{inherit}.{freq}.{region}.{score}.{het_comp}
            filter_combinations:
              - conservative.dominant.dominant_freq.all_genes.coding.passthrough

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
def variant_filtration_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return VariantFiltrationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_annotation": lambda x: "VAR_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return VariantFiltrationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FilterQualityStepPart ------------------------------------------------------------------


def test_filter_quality_step_part_get_input_files(variant_filtration_workflow):
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
        "VAR_ANNOTATION/output/bwa.gatk3_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1/out/"
        "bwa.gatk3_hc.jannovar_annotate_vcf.P001-N1-DNA1-WGS1"
    )
    root = get_project_root()
    pedigree_dict = {
        "ped": root + "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"
    }
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=var_base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_quality", "run")(wildcards)
    assert actual == expected


def test_filter_quality_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterQualityStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{thresholds}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_quality", "run")
    assert actual == expected


def test_filter_quality_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterQualityStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered."
        r"{index_library,[^\.]+}.{thresholds,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{thresholds}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_quality", "run")
    assert actual == expected


def test_filter_quality_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterQualityStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_quality", "run", resource)()
        assert actual == expected, msg_error


# Tests for FilterInheritanceStepPart   ------------------------------------------------------------


def test_filter_inheritance_step_part_get_input_files(variant_filtration_workflow):
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
        "work/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative/out/"
        "bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_inheritance", "run")(wildcards)
    assert actual == expected


def test_filter_inheritance_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}/out/{mapper}.{caller}."
        r"jannovar_annotate_vcf.filtered.{index_library}.{thresholds}.{inheritance}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_inheritance", "run")
    assert actual == expected


def test_filter_inheritance_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf."
        r"filtered.{index_library}.{thresholds}.{inheritance}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_inheritance", "run")
    assert actual == expected


def test_filter_inheritance_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterInheritanceStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_inheritance", "run", resource)()
        assert actual == expected, msg_error


# Tests for FilterFrequencyStepPart   --------------------------------------------------------------


def test_filter_frequency_step_part_get_input_files(variant_filtration_workflow):
    """Tests FilterFrequencyStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
            "inheritance": "dominant",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative.dominant/"
        "out/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative.dominant"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_frequency", "run")(wildcards)
    assert actual == expected


def test_filter_frequency_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterFrequencyStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{thresholds}."
        r"{inheritance}.{frequency}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_frequency", "run")
    assert actual == expected


def test_filter_frequency_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterFrequencyStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}/out/{mapper}.{caller}."
        r"jannovar_annotate_vcf.filtered.{index_library}.{thresholds}.{inheritance}.{frequency}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_frequency", "run")
    assert actual == expected


def test_filter_frequency_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterFrequencyStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_frequency", "run", resource)()
        assert actual == expected, msg_error


# Tests for FilterRegionsStepPart   ----------------------------------------------------------------


def test_filter_regions_step_part_get_input_files(variant_filtration_workflow):
    """Tests FilterRegionsStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
            "inheritance": "dominant",
            "frequency": "af_dominant",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative.dominant."
        "af_dominant/out/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1."
        "conservative.dominant.af_dominant"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_regions", "run")(wildcards)
    assert actual == expected


def test_filter_regions_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterRegionsStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{thresholds}."
        r"{inheritance}.{frequency}.{regions}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_regions", "run")
    assert actual == expected


def test_filter_regions_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterRegionsStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}/out/"
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{thresholds}."
        r"{inheritance}.{frequency}.{regions}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_regions", "run")
    assert actual == expected


def test_filter_regions_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterRegionsStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_regions", "run", resource)()
        assert actual == expected, msg_error


# Tests for FilterScoresStepPart   -----------------------------------------------------------------


def test_filter_scores_step_part_get_input_files(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
            "inheritance": "dominant",
            "frequency": "af_dominant",
            "regions": "all_genes",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative."
        "dominant.af_dominant.all_genes/out/bwa.gatk3_hc.jannovar_annotate_vcf.filtered."
        "P001-N1-DNA1-WGS1.conservative.dominant.af_dominant.all_genes"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_scores", "run")(wildcards)
    assert actual == expected


def test_filter_scores_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}."
        r"{thresholds}.{inheritance}.{frequency}.{regions}.{scores}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_scores", "run")
    assert actual == expected


def test_filter_scores_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}."
        r"{thresholds}.{inheritance}.{frequency}.{regions}.{scores}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_scores", "run")
    assert actual == expected


def test_filter_scores_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_scores", "run", resource)()
        assert actual == expected, msg_error


# Tests for FilterScoresStepPart   -----------------------------------------------------------------


def test_filter_het_comp_step_part_get_input_files(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "gatk3_hc",
            "index_library": "P001-N1-DNA1-WGS1",
            "thresholds": "conservative",
            "inheritance": "dominant",
            "frequency": "af_dominant",
            "regions": "all_genes",
            "scores": "coding",
        }
    )
    # Define expected
    base_name = (
        "work/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative."
        "dominant.af_dominant.all_genes.coding/out/bwa.gatk3_hc.jannovar_annotate_vcf."
        "filtered.P001-N1-DNA1-WGS1.conservative.dominant.af_dominant.all_genes.coding"
    )
    pedigree_dict = {"ped": "/work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    var_filtration_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    expected = {**pedigree_dict, **var_filtration_dict}
    # Get actual
    actual = variant_filtration_workflow.get_input_files("filter_het_comp", "run")(wildcards)
    assert actual == expected


def test_filter_het_comp_step_part_get_output_files(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}.{het_comp,[^\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf."
        r"filtered.{index_library}.{thresholds}.{inheritance}.{frequency}.{regions}.{scores}."
        r"{het_comp}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_filtration_workflow.get_output_files("filter_het_comp", "run")
    assert actual == expected


def test_filter_het_comp_step_part_get_log_file(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}.{het_comp,[^\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf."
        r"filtered.{index_library}.{thresholds}.{inheritance}.{frequency}.{regions}.{scores}."
        r"{het_comp}.log"
    )
    actual = variant_filtration_workflow.get_log_file("filter_het_comp", "run")
    assert actual == expected


def test_filter_het_comp_step_part_get_resource_usage(variant_filtration_workflow):
    """Tests FilterScoresStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_filtration_workflow.get_resource("filter_het_comp", "run", resource)()
        assert actual == expected, msg_error


# Tests for VariantFiltrationWorkflow  -------------------------------------------------------------


def test_variant_filtration_workflow(variant_filtration_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "filter_frequency",
        "filter_het_comp",
        "filter_inheritance",
        "filter_quality",
        "filter_regions",
        "filter_scores",
        "link_out",
        "write_pedigree",
    ]
    actual = list(sorted(variant_filtration_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P00{i}-N1-DNA1-WGS1.conservative."
        "dominant.dominant_freq.all_genes.coding.passthrough/out/"
        "bwa.gatk3_hc.jannovar_annotate_vcf.filtered.P00{i}-N1-DNA1-WGS1.conservative."
        "dominant.dominant_freq.all_genes.coding.passthrough.{ext}"
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
    actual = sorted(variant_filtration_workflow.get_result_files())
    assert actual == expected
