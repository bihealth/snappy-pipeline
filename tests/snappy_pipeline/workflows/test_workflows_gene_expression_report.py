# -*- coding: utf-8 -*-
"""Tests for the gene_expression_report workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.gene_expression_report import GeneExpressionReportWorkflow

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
              rna: ['star']
            star:
              path_index: /path/to/star/index

          gene_expression_quantification:
            tools: ['strandedness']

          gene_expression_report:
            path_gene_expression_quantification: /work

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
def gene_expression_report_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return GeneExpressionReportWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "gene_expression_quantification": lambda x: "GENE_EXPRESSION_QUANTIFICATION/" + x,
    }
    # Construct the workflow object
    return GeneExpressionReportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Global tests -------------------------------------------------------------------------------------


def test_all_steps_get_resource_usage(gene_expression_report_workflow):
    """Tests get_resource() for all steps"""
    # All available steps in workflow - all contain a single action, 'run'
    steps = (
        "aggregate_counts",
        "compute_ranks",
        "compute_signatures",
        "plot_expression_distribution",
    )
    # Define expected: default defined in workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for step in steps:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in step '{step}'."
            actual = gene_expression_report_workflow.get_resource(step, "run", resource)
            assert actual == expected, msg_error


# Tests for GeneExpressionReportAggreateFeaturecounts ----------------------------------------------


def test_gene_expression_rep_aggregate_feature_counts_step_part_get_input_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportAggreateFeaturecounts.get_input_files()"""
    base_out = (
        "GENE_EXPRESSION_QUANTIFICATION/output/star.featurecounts.P00{i}-T{t}-RNA1-mRNA_seq1/out/"
        "star.featurecounts.P00{i}-T{t}-RNA1-mRNA_seq1.tsv"
    )
    expected = [base_out.format(i=1, t=1), base_out.format(i=2, t=2)]
    actual = gene_expression_report_workflow.get_input_files("aggregate_counts", "run")
    assert actual == expected


def test_gene_expression_rep_aggregate_feature_counts_step_part_get_output_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportAggreateFeaturecounts.get_output_files()"""
    expected = {"tsv": "work/gene_exp.tsv"}
    actual = gene_expression_report_workflow.get_output_files("aggregate_counts", "run")
    assert actual == expected


def test_gene_expression_rep_aggregate_feature_counts_step_part_get_log_file(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportAggreateFeaturecounts.get_log_file()"""
    expected = (
        "work/{mapper}.aggregate_counts.{ngs_library}/log/"
        "snakemake.gene_expression_quantification.log"
    )
    actual = gene_expression_report_workflow.get_log_file("aggregate_counts", "run")
    assert actual == expected


# Tests for GeneExpressionReportRankExpression -----------------------------------------------------


def test_gene_expression_rep_compute_ranks_step_part_get_input_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportRankExpression.get_input_files()"""
    expected = {"tsv": "work/gene_exp.tsv"}
    actual = gene_expression_report_workflow.get_input_files("compute_ranks", "run")
    assert actual == expected


def test_gene_expression_rep_compute_ranks_step_part_get_output_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportRankExpression.get_output_files()"""
    expected = {"tsv": "work/{mapper}.{tool}.{ngs_library}/out/{mapper}.{tool}.{ngs_library}.tsv"}
    actual = gene_expression_report_workflow.get_output_files("compute_ranks", "run")
    assert actual == expected


def test_gene_expression_rep_compute_ranks_step_part_get_log_file(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportRankExpression.get_log_file()"""
    expected = (
        "work/{mapper}.compute_ranks.{ngs_library}/log/snakemake.gene_expression_quantification.log"
    )
    actual = gene_expression_report_workflow.get_log_file("compute_ranks", "run")
    assert actual == expected


# Tests for GeneExpressionReportComputeSignatures --------------------------------------------------


def test_gene_expression_rep_compute_signatures_step_part_get_output_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportComputeSignatures.get_output_files()"""
    expected = {"pdf": ["work/{mapper}.{tool}.{ngs_library}/out/{mapper}.{tool}.{ngs_library}.pdf"]}
    actual = gene_expression_report_workflow.get_output_files("compute_signatures", "run")
    assert actual == expected


def test_gene_expression_rep_compute_signatures_step_part_get_log_file(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportComputeSignatures.get_log_file()"""
    expected = (
        "work/{mapper}.compute_signatures.{ngs_library}/log/"
        "snakemake.gene_expression_quantification.log"
    )
    actual = gene_expression_report_workflow.get_log_file("compute_signatures", "run")
    assert actual == expected


# Tests for GeneExpressionReportPlotGeneDistribution -----------------------------------------------


def test_gene_expression_rep_plot_expression_distribution_step_part_get_output_files(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportPlotGeneDistribution.get_output_files()"""
    expected = {
        "pdf": ["work/{mapper}.{tool}.{ngs_library}/out/{mapper}.{tool}.{ngs_library}.genes.pdf"]
    }
    actual = gene_expression_report_workflow.get_output_files("plot_expression_distribution", "run")
    assert actual == expected


def test_gene_expression_rep_plot_expression_distribution_step_part_get_log_file(
    gene_expression_report_workflow,
):
    """Tests GeneExpressionReportPlotGeneDistribution.get_log_file()"""
    expected = (
        "work/{mapper}.plot_expression_distribution.{ngs_library}/log/"
        "snakemake.gene_expression_quantification.log"
    )
    actual = gene_expression_report_workflow.get_log_file("plot_expression_distribution", "run")
    assert actual == expected


# Tests for GeneExpressionReportWorkflow  ----------------------------------------------------------


def test_gene_expression_report_workflow(gene_expression_report_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "aggregate_counts",
        "compute_ranks",
        "compute_signatures",
        "link_out",
        "plot_expression_distribution",
    ]
    actual = list(sorted(gene_expression_report_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    base_out = (
        "output/star.featurecounts.P00{i}-T{t}-RNA1-mRNA_seq1/out/"
        "star.featurecounts.P00{i}-T{t}-RNA1-mRNA_seq1.{ext}"
    )
    expected = [
        base_out.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 2))
        for ext in (
            "genes.pdf",
            "pdf",
            "tsv",
        )
    ]
    expected = set(expected)
    actual = set(gene_expression_report_workflow.get_result_files())
    assert actual == expected
