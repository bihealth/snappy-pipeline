# -*- coding: utf-8 -*-
"""Tests for the helper_gcnv_model_wgs workflow module code"""

import textwrap
import unittest.mock as mock

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.helper_gcnv_model_wgs import HelperBuildWgsGcnvModelWorkflow

from .conftest import patch_module_fs


class MockRule:
    """Mocks Snakemake.Rule"""

    output = ["work/bwa.gcnv_scatter_intervals.default/out/bwa.gcnv_scatter_intervals.default/"]


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

          gcnv:
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
    """Return HelperBuildGcnvModelWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return HelperBuildWgsGcnvModelWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Test for BuildGcnvModelStepPart ------------------------------------------------------------------


def test_gcnv_get_cnv_model_result_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart.get_cnv_model_result_files()"""
    # Define expected
    interval_file = (
        "work/bwa.gcnv_filter_intervals.default/out/bwa.gcnv_filter_intervals.default.interval_list"
    )
    ploidy_file = "work/bwa.gcnv_contig_ploidy.default/out/bwa.gcnv_contig_ploidy.default/.done"
    expected = sorted([interval_file, ploidy_file])
    # Get actual
    actual = helper_gcnv_model_workflow.substep_getattr("gcnv", "get_cnv_model_result_files")(None)
    actual = sorted(actual)
    assert actual == expected


def test_gcnv_get_resource(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart.get_resource()"""
    high_resource_action_list = (
        "call_cnvs_cohort_mode",
        "call_cnvs_case_mode",
        "post_germline_calls",
        "post_germline_calls_cohort_mode",
        "post_germline_calls_case_mode",
    )
    actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "scatter_intervals",
        "call_cnvs_cohort_mode",
        "post_germline_calls",
        "post_germline_calls_cohort_mode",
        "merge_cohort_vcfs",
    )
    expected_low = {
        "threads": 1,
        "time": "04:00:00",
        "memory": "7680M",
    }
    expected_high = {
        "threads": 16,
        "time": "2-00:00:00",
        "memory": "46080M",
    }
    for action in actions:
        for resource in expected_low.keys():
            if action == "filter_intervals" and resource == "memory":
                actual = helper_gcnv_model_workflow.get_resource("gcnv", action, resource)(
                    None, attempt=1
                )
                assert actual == "20480M"
                actual = helper_gcnv_model_workflow.get_resource("gcnv", action, resource)(
                    None, attempt=2
                )
                assert actual == "24576M"
                actual = helper_gcnv_model_workflow.get_resource("gcnv", action, resource)(
                    None, attempt=3
                )
                assert actual == "28672M"
            else:
                if action in high_resource_action_list:
                    actual = helper_gcnv_model_workflow.get_resource("gcnv", action, resource)
                    assert actual == expected_high.get(resource)
                else:
                    actual = helper_gcnv_model_workflow.get_resource("gcnv", action, resource)
                    assert actual == expected_low.get(resource)


def test_gcnv_call_cnvs_cohort_mode_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs_cohort_mode()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Define expected
    interval_list_shard_out = (
        "work/{mapper}.gcnv_scatter_intervals.default/out/"
        "{mapper}.gcnv_scatter_intervals.default/temp_{shard}/scattered.interval_list"
    )
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = "work/bwa.gcnv_contig_ploidy.default/out/bwa.gcnv_contig_ploidy.default/.done"
    intervals_out = "work/gcnv_annotate_gc.default/out/gcnv_annotate_gc.default.tsv"
    expected = {
        "interval_list_shard": interval_list_shard_out,
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
        "intervals": intervals_out,
    }
    # Get actual
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "call_cnvs_cohort_mode")(wildcards)
    assert actual == expected


def test_gcnv_get_input_files_post_germline_calls_cohort_mode(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs_cohort_mode()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "calls": [],
        "ploidy": "work/bwa.gcnv_contig_ploidy.default/out/bwa.gcnv_contig_ploidy.default/.done",
    }
    with mock.patch("snakemake.checkpoints") as patched_checkpoints:
        # Patch checkpoint
        patched_checkpoints.build_gcnv_model_scatter_intervals = MockCheckpoint()
        # Get actual
        actual = helper_gcnv_model_workflow.get_input_files(
            "gcnv", "post_germline_calls_cohort_mode"
        )(wildcards, patched_checkpoints)
        assert actual == expected


# Test for HelperBuildWgsGcnvModelWorkflow  --------------------------------------------------------


def test_helper_gcnv_model_workflow(helper_gcnv_model_workflow):
    """Tests HelperBuildWgsGcnvModelWorkflow.get_result_files()"""
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        pattern_out.format(i=i, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for ext in (
            "interval.vcf.gz",
            "ratio.tsv",
            "vcf.gz",
        )
    ]
    actual = sorted(helper_gcnv_model_workflow.get_result_files())
    assert actual == expected
