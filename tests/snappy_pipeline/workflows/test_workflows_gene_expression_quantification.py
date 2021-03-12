# -*- coding: utf-8 -*-
"""Tests for the hla_typing workflow module code"""

import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.gene_expression_quantification import (
    GeneExpressionQuantificationWorkflow,
)

from .conftest import patch_module_fs

__author__ = "Clemens Messerschmidt"


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
              rna: ['star']
            star:
              path_index: /path/to/star/index

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
def gene_expression_quantification_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return GeneExpressionQuantificationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return GeneExpressionQuantificationWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FeatureCountsStepPart ------------------------------------------------------------------


def test_featurecounts_step_part_get_input_files(gene_expression_quantification_workflow):
    # Define expected
    ngs_mapping_base_out = "NGS_MAPPING/output/star.P001-T1-RNA1-mRNA_seq1/out/"
    expected = {
        "bai": ngs_mapping_base_out + "star.P001-T1-RNA1-mRNA_seq1.bam.bai",
        "bam": ngs_mapping_base_out + "star.P001-T1-RNA1-mRNA_seq1.bam",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-RNA1-mRNA_seq1", "mapper": "star"})
    actual = gene_expression_quantification_workflow.get_input_files("featurecounts", "run")(
        wildcards
    )

    assert actual == expected


def test_featurecounts_step_part_get_output_files(gene_expression_quantification_workflow):
    # Define expected
    base_out = "work/{mapper}.featurecounts.{library_name}/out/"
    expected = {
        "tsv": base_out + "{mapper}.featurecounts.{library_name}.tsv",
        "tsv_md5": base_out + "{mapper}.featurecounts.{library_name}.tsv.md5",
        "summary": base_out + "{mapper}.featurecounts.{library_name}.tsv.summary",
        "summary_md5": base_out + "{mapper}.featurecounts.{library_name}.tsv.summary.md5",
    }

    # Get actual
    actual = gene_expression_quantification_workflow.get_output_files("featurecounts", "run")

    assert actual == expected


def test_featurecounts_step_part_get_log_file(gene_expression_quantification_workflow):
    """Tests if `get_log_file` provides the correct path."""
    expected = (
        "work/{mapper}.featurecounts.{library_name}/log/{mapper}.featurecounts.{library_name}.log"
    )
    actual = gene_expression_quantification_workflow.get_log_file("featurecounts", "run").get("log")
    assert actual == expected


# Tests for GeneExpressionQuantificationWorkflow ---------------------------------------------------


def test_gene_expression_quantification_workflow_substeps(gene_expression_quantification_workflow):
    """Tests simple functionality of the workflow: checks if sub steps are created,
    i.e., the tools associated with gene expression quantification."""
    # Check created sub steps
    expected = [
        "duplication",
        "dupradar",
        "featurecounts",
        "link_in",
        "link_out",
        "rnaseqc",
        "salmon",
        "stats",
        "strandedness",
    ]
    actual = list(sorted(gene_expression_quantification_workflow.sub_steps.keys()))
    assert actual == expected


def test_gene_expression_quantification_workflow_files(gene_expression_quantification_workflow):
    """Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: duplication, dupradar, featurecounts,
    link_in, link_out, rnaseqc, salmon, stats, strandedness."""
    # Check result file construction
    base_name_out = (
        "output/{step}.P00{i}-T{i}-RNA1-mRNA_seq1/out/{step}.P00{i}-T{i}-RNA1-mRNA_seq1{ext}"
    )
    # Add expected Salmon out files
    expected = [
        base_name_out.format(step="salmon", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".gene.sf",
            ".gene.sf.md5",
            ".transcript.sf",
            ".transcript.sf.md5",
        )
    ]
    # Add expected star duplication out files
    expected += [
        base_name_out.format(step="star.duplication", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".pos.DupRate.xls",
            ".pos.DupRate.xls.md5",
            ".seq.DupRate.xls",
            ".seq.DupRate.xls.md5",
        )
    ]
    # Add expected star dupradar out files
    expected += [
        base_name_out.format(step="star.dupradar", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".dupradar.tsv",
            ".dupradar.tsv.md5",
        )
    ]
    # Add expected star featurecounts out files
    expected += [
        base_name_out.format(step="star.featurecounts", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".tsv",
            ".tsv.md5",
            ".tsv.summary",
            ".tsv.summary.md5",
        )
    ]
    # Add expected star rnaseqc out files
    expected += [
        base_name_out.format(step="star.rnaseqc", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".{s}{e}".format(s=s, e=e)
            for s in (
                "gapLengthHist_high",
                "gapLengthHist_low",
                "gapLengthHist_medium",
                "meanCoverage_high",
                "meanCoverage_low",
                "meanCoverage_medium",
                "meanCoverageNorm_high",
                "meanCoverageNorm_low",
                "meanCoverageNorm_medium",
            )
            for e in (".txt", ".txt.md5")
        )
    ]
    # Add expected star rnaseqc out files - cont.
    expected += [
        base_name_out.format(step="star.rnaseqc", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".metrics.tsv",
            ".metrics.tsv.md5",
        )
    ]
    # Add expected star stats out files
    expected += [
        base_name_out.format(step="star.stats", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".read_alignment_report.tsv",
            ".read_alignment_report.tsv.md5",
        )
    ]
    # Add expected star strandedness out files
    expected += [
        base_name_out.format(step="star.strandedness", i=i, ext=ext)
        for i in (1, 2)  # only for indices
        for ext in (
            ".decision",
            ".decision.md5",
            ".tsv",
            ".tsv.md5",
        )
    ]
    expected = set(expected)

    # Get actual
    actual = set(gene_expression_quantification_workflow.get_result_files())

    assert actual == expected
