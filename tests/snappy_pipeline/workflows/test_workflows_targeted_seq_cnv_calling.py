# -*- coding: utf-8 -*-
"""Tests for the targeted_seq_cnv_calling workflow module code"""

from copy import deepcopy
import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.targeted_seq_cnv_calling import TargetedSeqCnvCallingWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

# List of valid actions - gCNV
GCNV_ACTIONS = [
    "preprocess_intervals",
    "annotate_gc",
    "filter_intervals",
    "scatter_intervals",
    "coverage",
    "contig_ploidy",
    "call_cnvs",
    "post_germline_calls",
    "merge_cohort_vcfs",
    "extract_ped",
]

# List of valid actions - XHMM
XHMM_ACTIONS = [
    "coverage",
    "merge_cov",
    "ref_stats",
    "filter_center",
    "pca",
    "normalize",
    "zscore_center",
    "refilter",
    "discover",
    "genotype",
    "extract_ped",
]


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (somatic) configuration"""
    return yaml.round_trip_load(
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

          targeted_seq_cnv_calling:
            tools:
              - xhmm
              - gcnv
            xhmm:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
            gcnv:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
                  path_uniquely_mapable_bed: /path/to/uniquely/mappable/variable/GRCh37/file.bed.gz

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
def targeted_seq_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return TargetedSeqCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return TargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Global tests -------------------------------------------------------------------------------------


def test_call_assertion(targeted_seq_cnv_calling_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(Exception):
        targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "_undefined_action_")


def test_target_seq_cnv_calling_workflow_files(targeted_seq_cnv_calling_workflow):
    """Tests TargetedSeqCnvCallingWorkflow::get_result_files()

    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: gCNV and XHMM.
    """
    # Define expected
    pattern_out = "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/out/bwa.{tool}.P00{i}-N1-DNA1-WGS1.{ext}"
    expected = [
        pattern_out.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv", "xhmm")
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
        )
    ]
    expected = sorted(expected)
    # Get actual
    actual = sorted(targeted_seq_cnv_calling_workflow.get_result_files())
    assert actual == expected


# Global GcnvStepPart Tests ------------------------------------------------------------------------


def test_gcnv_update_cluster_config(targeted_seq_cnv_calling_workflow, dummy_cluster_config):
    """Tests GcnvStepPart::update_cluster_config for all actions"""
    # Define expected
    expected = {"mem", "time", "ntasks"}
    # Iterate over
    for action in GCNV_ACTIONS:
        key = "targeted_seq_cnv_calling_gcnv_{0}".format(action)
        actual = set(dummy_cluster_config[key].keys())
        assert actual == expected


def test_gcnv_get_params(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_params for all actions"""
    for action in GCNV_ACTIONS:
        if action == "coverage":
            targeted_seq_cnv_calling_workflow.get_params("gcnv", action)
        else:
            with pytest.raises(Exception):
                targeted_seq_cnv_calling_workflow.get_params("gcnv", action)


# Tests for GcnvStepPart (preprocess_intervals) ----------------------------------------------------


def test_gcnv_preprocess_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_preprocess_intervals()"""
    # Define expected - empty dictionary for all
    expected = {}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "preprocess_intervals")(None)
    assert actual == expected


def test_gcnv_preprocess_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_preprocess_intervals()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.{library_kit}/out/"
        "gcnv_preprocess_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "preprocess_intervals")
    assert actual == expected


def test_gcnv_target_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'preprocess_intervals' step"""
    # Define expected
    expected = (
        "work/gcnv_preprocess_intervals.{library_kit}/log/"
        "gcnv_preprocess_intervals.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "preprocess_intervals")
    assert actual == expected


# Tests for GcnvStepPart (annotate_gc) -------------------------------------------------------------


def test_gcnv_annotate_gc_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_annotate_gc()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual - Note: library kit defined in conftest: germline_sheet_tsv
    wildcards = Wildcards(fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "annotate_gc")(wildcards)
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_annotate_gc()"""
    # Define expected
    expected = {"tsv": "work/gcnv_annotate_gc.{library_kit}/out/gcnv_annotate_gc.{library_kit}.tsv"}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "annotate_gc")
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'annotate_gc' step"""
    # Define expected
    expected = "work/gcnv_annotate_gc.{library_kit}/log/gcnv_annotate_gc.{library_kit}.log"
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "annotate_gc")
    assert actual == expected


# Tests for GcnvStepPart (filter_intervals) --------------------------------------------------------


def test_gcnv_filter_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_filter_intervals()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "filter_intervals")(
        wildcards
    )
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_filter_intervals()"""
    # Define expected
    interval_list_out = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": interval_list_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "filter_intervals")
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'filter_intervals' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/log/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "filter_intervals")
    assert actual == expected


# Tests for GcnvStepPart (scatter_intervals) -------------------------------------------------------


def test_gcnv_scatter_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_scatter_intervals()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "scatter_intervals")(
        wildcards
    )
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_scatter_intervals()"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "scatter_intervals")
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'scatter_intervals' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/log/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "scatter_intervals")
    assert actual == expected


# Tests for GcnvStepPart (coverage) ----------------------------------------------------------------


def test_gcnv_coverage_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_coverage()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "coverage")(wildcards)
    assert actual == expected


def test_gcnv_coverage_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_coverage()"""
    # Define expected
    tsv_out = (
        "work/{mapper}.gcnv_coverage.{library_name}/out/{mapper}.gcnv_coverage.{library_name}.tsv"
    )
    expected = {"tsv": tsv_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "coverage")
    assert actual == expected


def test_gcnv_coverage_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'coverage' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_coverage.{library_name}/log/{mapper}.gcnv_coverage.{library_name}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "coverage")
    assert actual == expected


# Tests for GcnvStepPart (contig_ploidy) -----------------------------------------------------------


def test_gcnv_contig_ploidy_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_contig_ploidy()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_contig_ploidy()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "contig_ploidy")
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'contig_ploidy' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/log/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "contig_ploidy")
    assert actual == expected


# Tests for GcnvStepPart (call_cnvs) ---------------------------------------------------------------


def test_gcnv_call_cnvs_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_call_cnvs()"""
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
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "call_cnvs")(wildcards)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_call_cnvs()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "call_cnvs")
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'call_cnvs' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "call_cnvs")
    assert actual == expected


# Tests for GcnvStepPart (merge_cohort_vcfs) -----------------------------------------------------


def test_gcnv_merge_cohort_vcfs_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_merge_cohort_vcfs()"""
    # Define expected
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.vcf.gz"
    )
    expected = [pattern_out.format(i=i) for i in range(1, 7)]  # P001 - P006
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "merge_cohort_vcfs")(
        wildcards
    )
    assert actual == expected


def test_gcnv_merge_cohort_vcfs_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_merge_cohort_vcfs()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.gcnv_merge_cohort_vcfs.{library_kit}/out/"
        "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "merge_cohort_vcfs")
    assert actual == expected


def test_gcnv_merge_cohort_vcfs_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'merge_cohort_vcfs' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_merge_cohort_vcfs.{library_kit}/log/"
        "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "merge_cohort_vcfs")
    assert actual == expected


# Tests for GcnvStepPart (extract_ped) -------------------------------------------------------------


def test_gcnv_extract_ped_vcfs_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_input_files_extract_ped()"""
    # Define expected
    pattern_out = (
        "work/bwa.gcnv_merge_cohort_vcfs.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_merge_cohort_vcfs.Agilent_SureSelect_Human_All_Exon_V6"
    )
    expected = {"vcf": pattern_out + ".vcf.gz", "tbi": pattern_out + ".vcf.gz.tbi"}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "extract_ped")(wildcards)
    assert actual == expected


def test_gcnv_extract_ped_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::_get_output_files_extract_ped()"""
    # Define expected
    pattern_out = "work/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "extract_ped")
    assert actual == expected


def test_gcnv_extract_ped_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart::get_log_file for 'extract_ped' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_extract_ped.{library_name}/log/"
        "{mapper}.gcnv_extract_ped.{library_name}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "extract_ped")
    assert actual == expected


# Global XhmmStepPart Tests ------------------------------------------------------------------------


def test_xhmm_get_params(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart::get_params for all actions"""
    for action in XHMM_ACTIONS:
        if action == "coverage":
            targeted_seq_cnv_calling_workflow.get_params("xhmm", action)
        else:
            with pytest.raises(Exception):
                targeted_seq_cnv_calling_workflow.get_params("xhmm", action)
