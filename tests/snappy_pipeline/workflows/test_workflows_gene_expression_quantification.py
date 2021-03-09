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
                rna:
                    - star
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


# Tests for FeatureCountsStepPart ----------------------------------------------------------------------


def test_featurecounts_step_part_get_input_files(gene_expression_quantification_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-RNA1-mRNA_seq1", "mapper": "star"})
    expected = {
        "bai": "NGS_MAPPING/output/star.P001-T1-RNA1-mRNA_seq1/out/star.P001-T1-RNA1-mRNA_seq1.bam.bai",
        "bam": "NGS_MAPPING/output/star.P001-T1-RNA1-mRNA_seq1/out/star.P001-T1-RNA1-mRNA_seq1.bam",
    }
    assert (
        gene_expression_quantification_workflow.get_input_files("featurecounts", "run")(wildcards)
        == expected
    )


def test_featurecounts_step_part_get_output_files(gene_expression_quantification_workflow):
    expected = {
        "tsv": "work/{mapper}.featurecounts.{library_name}/out/{mapper}.featurecounts.{library_name}.tsv",
        "tsv_md5": "work/{mapper}.featurecounts.{library_name}/out/{mapper}.featurecounts.{library_name}.tsv.md5",
        "summary": "work/{mapper}.featurecounts.{library_name}/out/{mapper}.featurecounts.{library_name}.tsv.summary",
        "summary_md5": "work/{mapper}.featurecounts.{library_name}/out/{mapper}.featurecounts.{library_name}.tsv.summary.md5",
    }
    assert (
        gene_expression_quantification_workflow.get_output_files("featurecounts", "run") == expected
    )


def test_featurecounts_step_part_get_log_file(gene_expression_quantification_workflow):
    """Tests if `get_log_file` provides the correct path."""
    expected = (
        "work/{mapper}.featurecounts.{library_name}/log/{mapper}.featurecounts.{library_name}.log"
    )
    actual = gene_expression_quantification_workflow.get_log_file("featurecounts", "run").get("log")
    print(actual)
    print(expected)
    assert actual == expected


# Tests for GeneExpressionQuantificationWorkflow ---------------------------------------------------------------------


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
    expected = [
        "output/salmon.P001-T1-RNA1-mRNA_seq1/out/salmon.P001-T1-RNA1-mRNA_seq1.gene.sf",
        "output/salmon.P001-T1-RNA1-mRNA_seq1/out/salmon.P001-T1-RNA1-mRNA_seq1.gene.sf.md5",
        "output/salmon.P001-T1-RNA1-mRNA_seq1/out/salmon.P001-T1-RNA1-mRNA_seq1.transcript.sf",
        "output/salmon.P001-T1-RNA1-mRNA_seq1/out/salmon.P001-T1-RNA1-mRNA_seq1.transcript.sf.md5",
        "output/salmon.P002-T2-RNA1-mRNA_seq1/out/salmon.P002-T2-RNA1-mRNA_seq1.gene.sf",
        "output/salmon.P002-T2-RNA1-mRNA_seq1/out/salmon.P002-T2-RNA1-mRNA_seq1.gene.sf.md5",
        "output/salmon.P002-T2-RNA1-mRNA_seq1/out/salmon.P002-T2-RNA1-mRNA_seq1.transcript.sf",
        "output/salmon.P002-T2-RNA1-mRNA_seq1/out/salmon.P002-T2-RNA1-mRNA_seq1.transcript.sf.md5",
        "output/star.duplication.P001-T1-RNA1-mRNA_seq1/out/star.duplication.P001-T1-RNA1-mRNA_seq1.pos.DupRate.xls",
        "output/star.duplication.P001-T1-RNA1-mRNA_seq1/out/star.duplication.P001-T1-RNA1-mRNA_seq1.pos.DupRate.xls.md5",
        "output/star.duplication.P001-T1-RNA1-mRNA_seq1/out/star.duplication.P001-T1-RNA1-mRNA_seq1.seq.DupRate.xls",
        "output/star.duplication.P001-T1-RNA1-mRNA_seq1/out/star.duplication.P001-T1-RNA1-mRNA_seq1.seq.DupRate.xls.md5",
        "output/star.duplication.P002-T2-RNA1-mRNA_seq1/out/star.duplication.P002-T2-RNA1-mRNA_seq1.pos.DupRate.xls",
        "output/star.duplication.P002-T2-RNA1-mRNA_seq1/out/star.duplication.P002-T2-RNA1-mRNA_seq1.pos.DupRate.xls.md5",
        "output/star.duplication.P002-T2-RNA1-mRNA_seq1/out/star.duplication.P002-T2-RNA1-mRNA_seq1.seq.DupRate.xls",
        "output/star.duplication.P002-T2-RNA1-mRNA_seq1/out/star.duplication.P002-T2-RNA1-mRNA_seq1.seq.DupRate.xls.md5",
        "output/star.dupradar.P001-T1-RNA1-mRNA_seq1/out/star.dupradar.P001-T1-RNA1-mRNA_seq1.dupradar.tsv",
        "output/star.dupradar.P001-T1-RNA1-mRNA_seq1/out/star.dupradar.P001-T1-RNA1-mRNA_seq1.dupradar.tsv.md5",
        "output/star.dupradar.P002-T2-RNA1-mRNA_seq1/out/star.dupradar.P002-T2-RNA1-mRNA_seq1.dupradar.tsv",
        "output/star.dupradar.P002-T2-RNA1-mRNA_seq1/out/star.dupradar.P002-T2-RNA1-mRNA_seq1.dupradar.tsv.md5",
        "output/star.featurecounts.P001-T1-RNA1-mRNA_seq1/out/star.featurecounts.P001-T1-RNA1-mRNA_seq1.tsv",
        "output/star.featurecounts.P001-T1-RNA1-mRNA_seq1/out/star.featurecounts.P001-T1-RNA1-mRNA_seq1.tsv.md5",
        "output/star.featurecounts.P001-T1-RNA1-mRNA_seq1/out/star.featurecounts.P001-T1-RNA1-mRNA_seq1.tsv.summary",
        "output/star.featurecounts.P001-T1-RNA1-mRNA_seq1/out/star.featurecounts.P001-T1-RNA1-mRNA_seq1.tsv.summary.md5",
        "output/star.featurecounts.P002-T2-RNA1-mRNA_seq1/out/star.featurecounts.P002-T2-RNA1-mRNA_seq1.tsv",
        "output/star.featurecounts.P002-T2-RNA1-mRNA_seq1/out/star.featurecounts.P002-T2-RNA1-mRNA_seq1.tsv.md5",
        "output/star.featurecounts.P002-T2-RNA1-mRNA_seq1/out/star.featurecounts.P002-T2-RNA1-mRNA_seq1.tsv.summary",
        "output/star.featurecounts.P002-T2-RNA1-mRNA_seq1/out/star.featurecounts.P002-T2-RNA1-mRNA_seq1.tsv.summary.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_high.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_high.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_low.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_low.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_medium.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.gapLengthHist_medium.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_high.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_high.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_low.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_low.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_medium.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverage_medium.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_high.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_high.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_low.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_low.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_medium.txt",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.meanCoverageNorm_medium.txt.md5",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.metrics.tsv",
        "output/star.rnaseqc.P001-T1-RNA1-mRNA_seq1/out/star.rnaseqc.P001-T1-RNA1-mRNA_seq1.metrics.tsv.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_high.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_high.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_low.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_low.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_medium.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.gapLengthHist_medium.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_high.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_high.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_low.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_low.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_medium.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverage_medium.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_high.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_high.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_low.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_low.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_medium.txt",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.meanCoverageNorm_medium.txt.md5",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.metrics.tsv",
        "output/star.rnaseqc.P002-T2-RNA1-mRNA_seq1/out/star.rnaseqc.P002-T2-RNA1-mRNA_seq1.metrics.tsv.md5",
        "output/star.stats.P001-T1-RNA1-mRNA_seq1/out/star.stats.P001-T1-RNA1-mRNA_seq1.read_alignment_report.tsv",
        "output/star.stats.P001-T1-RNA1-mRNA_seq1/out/star.stats.P001-T1-RNA1-mRNA_seq1.read_alignment_report.tsv.md5",
        "output/star.stats.P002-T2-RNA1-mRNA_seq1/out/star.stats.P002-T2-RNA1-mRNA_seq1.read_alignment_report.tsv",
        "output/star.stats.P002-T2-RNA1-mRNA_seq1/out/star.stats.P002-T2-RNA1-mRNA_seq1.read_alignment_report.tsv.md5",
        "output/star.strandedness.P001-T1-RNA1-mRNA_seq1/out/star.strandedness.P001-T1-RNA1-mRNA_seq1.decision",
        "output/star.strandedness.P001-T1-RNA1-mRNA_seq1/out/star.strandedness.P001-T1-RNA1-mRNA_seq1.decision.md5",
        "output/star.strandedness.P001-T1-RNA1-mRNA_seq1/out/star.strandedness.P001-T1-RNA1-mRNA_seq1.tsv",
        "output/star.strandedness.P001-T1-RNA1-mRNA_seq1/out/star.strandedness.P001-T1-RNA1-mRNA_seq1.tsv.md5",
        "output/star.strandedness.P002-T2-RNA1-mRNA_seq1/out/star.strandedness.P002-T2-RNA1-mRNA_seq1.decision",
        "output/star.strandedness.P002-T2-RNA1-mRNA_seq1/out/star.strandedness.P002-T2-RNA1-mRNA_seq1.decision.md5",
        "output/star.strandedness.P002-T2-RNA1-mRNA_seq1/out/star.strandedness.P002-T2-RNA1-mRNA_seq1.tsv",
        "output/star.strandedness.P002-T2-RNA1-mRNA_seq1/out/star.strandedness.P002-T2-RNA1-mRNA_seq1.tsv.md5",
    ]
    expected = set(expected)
    actual = set(gene_expression_quantification_workflow.get_result_files())
    assert actual == expected
