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


# Test for BuildGcnvModelStepPart ------------------------------------------------------------------


def test_gcnv_call_cnvs_cohort_mode_step_part_get_input_files(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs_cohort_mode()"""
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
    actual = helper_gcnv_model_workflow.get_input_files("gcnv", "call_cnvs_cohort_mode")(wildcards)
    assert actual == expected


def test_gcnv_get_input_files_post_germline_calls_cohort_mode(helper_gcnv_model_workflow):
    """Tests BuildGcnvModelStepPart._get_input_files_call_cnvs_cohort_mode()"""
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
        actual = helper_gcnv_model_workflow.get_input_files(
            "gcnv", "post_germline_calls_cohort_mode"
        )(wildcards, patched_checkpoints)
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
        for i in (1, 4)  # only index: P001, P004
        for ext in (
            "interval.vcf.gz",
            "ratio.tsv",
            "vcf.gz",
        )
    ]
    actual = sorted(helper_gcnv_model_workflow.get_result_files())
    assert actual == expected
