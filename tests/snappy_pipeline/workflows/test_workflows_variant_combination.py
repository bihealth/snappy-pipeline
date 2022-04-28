# -*- coding: utf-8 -*-
"""Tests for the variant_checking workflow module code"""

import textwrap

import pytest
from ruamel import yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_combination import VariantCombinationWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
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
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk_hc

          variant_combination:
            combinations:
              # Note: `all` added just for test purposes, it is commented out in original code
              - name: all
                operation: vars_intersect
                left: wgs_sv_filtration:no_filter.all.wgs
                right: variant_filtration:no_filter.de_novo.freq_all.wgs.score_all.passthrough

              # CNV and a conserved small variant share a common *limb* TAD.
              - name: small_conserved_cnv_common_limb_tad
                operation: vars_share_interval
                args:
                  intervals_bed: /path/to/static_data/GRCh37/newlimb_tads.bed
                left: wgs_cnv_filtration:conservative.dominant.whole_genome
                right: "variant_filtration:conservative.dominant.recessive_freq.wgs.conserved.\
                        passthrough"

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
def variant_combination_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return VariantCombinationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_filtration": lambda x: "VAR_FILTRATION/" + x,
        "wgs_sv_filtration": lambda x: "WGS_VAR_FILTRATION/" + x,
        "wgs_cnv_filtration": lambda x: "WGS_CNV_FILTRATION/" + x,
    }
    # Construct the workflow object
    return VariantCombinationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VarsIntersectStepPart ------------------------------------------------------------------


def test_vars_intersect_step_part_get_input_files(variant_combination_workflow):
    """Tests VarsIntersectStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "index_library": "P001-N1-DNA1-WGS1",
            "combination": "all",
            "left_caller": "wgs_sv_filtration",
            "right_caller": "variant_filtration",
        }
    )
    # Define expected
    wgs_base_name = (
        "WGS_VAR_FILTRATION/output/"
        "bwa.wgs_sv_filtration.annotated.filtered.P001-N1-DNA1-WGS1.no_filter.all.wgs/out/"
        "bwa.wgs_sv_filtration.annotated.filtered.P001-N1-DNA1-WGS1.no_filter.all.wgs"
    )
    var_base_name = (
        "VAR_FILTRATION/output/"
        "bwa.variant_filtration.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.no_filter."
        "de_novo.freq_all.wgs.score_all.passthrough/out/bwa.variant_filtration."
        "jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.no_filter.de_novo.freq_all."
        "wgs.score_all.passthrough"
    )
    pedigree_dict = {"ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    _tmp_wgs_filtration_dict = get_expected_output_vcf_files_dict(base_out=wgs_base_name)
    wgs_filtration_dict = {}
    for key, value in _tmp_wgs_filtration_dict.items():
        new_key = "left_" + key
        wgs_filtration_dict[new_key] = value
    _tmp_var_filtration_dict = get_expected_output_vcf_files_dict(base_out=var_base_name)
    var_filtration_dict = {}
    for key, value in _tmp_var_filtration_dict.items():
        new_key = "right_" + key
        var_filtration_dict[new_key] = value
    expected = {**pedigree_dict, **wgs_filtration_dict, **var_filtration_dict}
    # Get actual
    actual = variant_combination_workflow.get_input_files("vars_intersect", "run")(wildcards)
    assert actual == expected


def test_vars_intersect_step_part_get_output_files(variant_combination_workflow):
    """Tests VarsIntersectStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.combined_variants.vars_intersect.{index_library}.{combination}."
        "{left_caller}.{right_caller}/out/{mapper}.combined_variants.vars_intersect."
        "{index_library}.{combination}.{left_caller}.{right_caller}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_combination_workflow.get_output_files("vars_intersect", "run")
    assert actual == expected


def test_vars_intersect_step_part_get_log_file(variant_combination_workflow):
    """Tests VarsIntersectStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.combined_variants.vars_intersect.{index_library}.{combination}."
        "{left_caller}.{right_caller}/out/{mapper}.combined_variants.vars_intersect."
        "{index_library}.{combination}.{left_caller}.{right_caller}.log"
    )
    actual = variant_combination_workflow.get_log_file("vars_intersect", "run")
    assert actual == expected


def test_vars_intersect_step_part_get_args(variant_combination_workflow):
    """Tests VarsIntersectStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"combination": "all"})
    expected = {}
    actual = variant_combination_workflow.get_args("vars_intersect", "run")(wildcards)
    assert actual == expected


def test_vars_intersect_step_part_get_resource_usage(variant_combination_workflow):
    """Tests VarsIntersectStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7577M", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_combination_workflow.get_resource("vars_intersect", "run", resource)
        assert actual == expected, msg_error


# Tests for VarsShareIntervalStepPart --------------------------------------------------------------


def test_vars_share_interval_step_part_get_input_files(variant_combination_workflow):
    """Tests VarsShareIntervalStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "index_library": "P001-N1-DNA1-WGS1",
            "combination": "small_conserved_cnv_common_limb_tad",
            "left_caller": "wgs_cnv_filtration",
            "right_caller": "variant_filtration",
        }
    )
    # Define expected
    wgs_base_name = (
        "WGS_CNV_FILTRATION/output/bwa.wgs_cnv_filtration.annotated.filtered.P001-N1-DNA1-WGS1."
        "conservative.dominant.whole_genome/out/bwa.wgs_cnv_filtration.annotated.filtered."
        "P001-N1-DNA1-WGS1.conservative.dominant.whole_genome"
    )
    var_base_name = (
        "VAR_FILTRATION/output/bwa.variant_filtration.jannovar_annotate_vcf.filtered."
        "P001-N1-DNA1-WGS1.conservative.dominant.recessive_freq.wgs.conserved.passthrough/out/"
        "bwa.variant_filtration.jannovar_annotate_vcf.filtered.P001-N1-DNA1-WGS1.conservative."
        "dominant.recessive_freq.wgs.conserved.passthrough"
    )
    pedigree_dict = {"ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    _tmp_wgs_filtration_dict = get_expected_output_vcf_files_dict(base_out=wgs_base_name)
    wgs_filtration_dict = {}
    for key, value in _tmp_wgs_filtration_dict.items():
        new_key = "left_" + key
        wgs_filtration_dict[new_key] = value
    _tmp_var_filtration_dict = get_expected_output_vcf_files_dict(base_out=var_base_name)
    var_filtration_dict = {}
    for key, value in _tmp_var_filtration_dict.items():
        new_key = "right_" + key
        var_filtration_dict[new_key] = value
    expected = {**pedigree_dict, **wgs_filtration_dict, **var_filtration_dict}
    # Get actual
    actual = variant_combination_workflow.get_input_files("vars_share_interval", "run")(wildcards)
    assert actual == expected


def test_vars_share_interval_step_part_get_output_files(variant_combination_workflow):
    """Tests VarsShareIntervalStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.combined_variants.vars_share_interval.{index_library}.{combination}."
        "{left_caller}.{right_caller}/out/{mapper}.combined_variants.vars_share_interval."
        "{index_library}.{combination}.{left_caller}.{right_caller}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_combination_workflow.get_output_files("vars_share_interval", "run")
    assert actual == expected


def test_vars_share_interval_step_part_get_log_file(variant_combination_workflow):
    """Tests VarsShareIntervalStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.combined_variants.vars_share_interval.{index_library}.{combination}."
        "{left_caller}.{right_caller}/out/{mapper}.combined_variants.vars_share_interval."
        "{index_library}.{combination}.{left_caller}.{right_caller}.log"
    )
    actual = variant_combination_workflow.get_log_file("vars_share_interval", "run")
    assert actual == expected


def test_vars_share_interval_step_part_get_args(variant_combination_workflow):
    """Tests VarsShareIntervalStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"combination": "small_conserved_cnv_common_limb_tad"})
    expected = {"intervals_bed": "/path/to/static_data/GRCh37/newlimb_tads.bed"}
    actual = variant_combination_workflow.get_args("vars_share_interval", "run")(wildcards)
    assert actual == expected


def test_vars_share_interval_step_part_get_resource_usage(variant_combination_workflow):
    """Tests VarsShareIntervalStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7577M", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_combination_workflow.get_resource("vars_share_interval", "run", resource)
        assert actual == expected, msg_error


# Tests for VariantCombinationWorkflow -------------------------------------------------------------


def test_variant_checking_workflow(variant_combination_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "vars_intersect", "vars_share_interval", "write_pedigree"]
    actual = list(sorted(variant_combination_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/bwa.combined_variants.{type_}."
        "P00{i}-N1-DNA1-WGS1.{combination}.{left}.{right}/out/"
        "bwa.combined_variants.{type_}.P00{i}-N1-DNA1-WGS1.{combination}.{left}.{right}.{ext}"
    )
    expected = [
        tpl.format(type_=type_, combination=combination, left=left, right=right, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
        )
        for type_, combination, left, right in (
            ("vars_intersect", "all", "delly2", "gatk_hc"),
            ("vars_share_interval", "small_conserved_cnv_common_limb_tad", "erds_sv2", "gatk_hc"),
        )
    ]
    expected = sorted(expected)
    actual = sorted(variant_combination_workflow.get_result_files())
    assert actual == expected
