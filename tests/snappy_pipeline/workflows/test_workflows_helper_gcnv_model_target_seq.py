# -*- coding: utf-8 -*-
"""Tests for the helper_gcnv_model_target_seq workflow module code"""

import textwrap
import unittest.mock as mock

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.helper_gcnv_model_target_seq import (
    HelperBuildTargetSeqGcnvModelWorkflow,
)

from .common import get_expected_gcnv_log_file
from .conftest import patch_module_fs


class MockRule:
    """Mocks Snakemake.Rule"""

    output = [
        (
            "work/bwa.gcnv_scatter_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
            "bwa.gcnv_scatter_intervals.Agilent_SureSelect_Human_All_Exon_V6/"
        )
    ]


class MockCheckpoint:
    """Mocks Snakemake.Checkpoint"""

    def get(self, **wildcards):
        _ = {**wildcards}
        return MockRule()


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
          helper_gcnv_model_target_seq:
            gcnv:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
              path_uniquely_mapable_bed: /path/to/map_track.bed  # REQUIRED

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
def helper_gcnv_model_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return HelperBuildTargetSeqGcnvModelWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return HelperBuildTargetSeqGcnvModelWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for BuildGcnvTargetSeqModelStepPart (preprocess_intervals) ---------------------------------


def test_gcnv_preprocess_intervals_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_preprocess_intervals()"""
    expected = {}
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "preprocess_intervals")(None)
    assert actual == expected


def test_gcnv_preprocess_intervals_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_preprocess_intervals()"""
    output_path = (
        "work/gcnv_preprocess_intervals.{library_kit}/out/"
        "gcnv_preprocess_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": output_path}
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "preprocess_intervals")
    assert actual == expected


def test_gcnv_target_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'preprocess_intervals' step"""
    expected = (
        "work/gcnv_preprocess_intervals.{library_kit}/log/"
        "gcnv_preprocess_intervals.{library_kit}.log"
    )
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "preprocess_intervals")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (coverage) ---------------------------------------------


def test_gcnv_coverage_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_coverage()"""
    # Define expected
    interval_list_out = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    bam_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1"
    expected = {
        "interval_list": interval_list_out,
        "bam": bam_out + ".bam",
        "bai": bam_out + ".bam.bai",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "coverage")(wildcards)
    assert actual == expected


def test_gcnv_coverage_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_coverage()"""
    tsv_out = (
        "work/{mapper}.gcnv_coverage.{library_name}/out/{mapper}.gcnv_coverage.{library_name}.tsv"
    )
    expected = {"tsv": tsv_out}
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "coverage")
    assert actual == expected


def test_gcnv_coverage_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'coverage' step"""
    expected = (
        "work/{mapper}.gcnv_coverage.{library_name}/log/{mapper}.gcnv_coverage.{library_name}.log"
    )
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "coverage")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (annotate_gc) ------------------------------------------


def test_gcnv_annotate_gc_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_annotate_gc()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual - Note: library kit defined in conftest: germline_sheet_tsv
    wildcards = Wildcards(fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6"})
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "annotate_gc")(wildcards)
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_annotate_gc()"""
    # Define expected
    expected = {"tsv": "work/gcnv_annotate_gc.{library_kit}/out/gcnv_annotate_gc.{library_kit}.tsv"}
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "annotate_gc")
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'annotate_gc' step"""
    # Define expected
    expected = "work/gcnv_annotate_gc.{library_kit}/log/gcnv_annotate_gc.{library_kit}.log"
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "annotate_gc")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (filter_intervals) -------------------------------------


def test_gcnv_filter_intervals_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_filter_intervals()"""
    # Define expected
    interval_list_out = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    tsv_out = (
        "work/gcnv_annotate_gc.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_annotate_gc.Agilent_SureSelect_Human_All_Exon_V6.tsv"
    )
    csv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    csv_list_out = [csv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {"interval_list": interval_list_out, "tsv": tsv_out, "covs": csv_list_out}
    # Get actual. Notes:
    # - library kit defined in conftest: `germline_sheet_tsv`
    # - mapper defined in `minimal_config`
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "filter_intervals")(wildcards)
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_filter_intervals()"""
    # Define expected
    interval_list_out = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": interval_list_out}
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "filter_intervals")
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'filter_intervals' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="filter_intervals")
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "filter_intervals")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (scatter_intervals) ------------------------------------


def test_gcnv_scatter_intervals_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_scatter_intervals()"""
    # Define expected
    output_path = (
        "work/bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual - Note: library kit defined in conftest: germline_sheet_tsv
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "scatter_intervals")(wildcards)
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_scatter_intervals()"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}"
    )
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "scatter_intervals")
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'scatter_intervals' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="scatter_intervals")
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "scatter_intervals")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (contig_ploidy) ----------------------------------------


def test_gcnv_contig_ploidy_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_input_files_contig_ploidy()"""
    # Define expected
    interval_list_out = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
    )
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {"interval_list": interval_list_out, "tsv": tsv_list_out}
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_contig_ploidy()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "contig_ploidy")
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'contig_ploidy' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="contig_ploidy")
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "contig_ploidy")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (call_cnvs) --------------------------------------------


def test_gcnv_call_cnvs_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs()"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    # Define expected
    interval_list_shard_out = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}/temp_{shard}/scattered.interval_list"
    )
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    intervals_out = "work/gcnv_annotate_gc.{library_kit}/out/gcnv_annotate_gc.{library_kit}.tsv"
    expected = {
        "interval_list_shard": interval_list_shard_out,
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
        "intervals": intervals_out,
    }
    # Get actual
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "call_cnvs")(wildcards)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_call_cnvs()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "call_cnvs")
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file for 'call_cnvs' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}.log"
    )
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "call_cnvs")
    assert actual == expected


# Tests for BuildGcnvTargetSeqModelStepPart (post_germline_calls) ----------------------------------


def test_gcnv_get_input_files_post_germline_calls(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "calls": [],
        "ploidy": (
            "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
            "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
        ),
    }
    with mock.patch("snakemake.checkpoints") as patched_checkpoints:
        # Patch checkpoint
        patched_checkpoints.build_gcnv_model_scatter_intervals = MockCheckpoint()
        # Get actual
        actual = helper_gcnv_model_workflow.get_input_files("gcnv", "post_germline_calls")(
            wildcards, patched_checkpoints
        )
        assert actual == expected


def test_gcnv_post_germline_calls_step_part_get_output_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart._get_output_files_post_germline_calls()"""
    # Define expected
    base_name = (
        "work/{mapper}.gcnv_post_germline_calls.{library_name}/out/"
        "{mapper}.gcnv_post_germline_calls.{library_name}"
    )
    expected = {
        "ratio_tsv": base_name + ".ratio.tsv",
        "itv_vcf": base_name + ".interval.vcf.gz",
        "seg_vcf": base_name + ".vcf.gz",
    }
    # Get actual
    actual = helper_gcnv_model_workflow.get_output_files("gcnv", "post_germline_calls")
    assert actual == expected


def test_gcnv_post_germline_calls_step_part_get_log_file(helper_gcnv_model_workflow):
    """Tests BuildGcnvTargetSeqModelStepPart.get_log_file() for 'post_germline_calls' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_post_germline_calls.{library_name}/log/"
        "{mapper}.gcnv_post_germline_calls.{library_name}.log"
    )
    # Get actual
    actual = helper_gcnv_model_workflow.get_log_file("gcnv", "post_germline_calls")
    assert actual == expected


# Test for HelperBuildTargetSeqGcnvModelWorkflow  --------------------------------------------------


def test_helper_gcnv_model_workflow(helper_gcnv_model_workflow):
    """Tests HelperBuildTargetSeqGcnvModelWorkflow.get_result_files()"""
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        pattern_out.format(i=i, ext=ext)
        for i in range(1, 7)
        for ext in (
            "interval.vcf.gz",
            "ratio.tsv",
            "vcf.gz",
        )
    ]
    actual = sorted(helper_gcnv_model_workflow.get_result_files())
    assert actual == expected
