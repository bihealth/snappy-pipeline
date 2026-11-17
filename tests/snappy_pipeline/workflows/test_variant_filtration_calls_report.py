# -*- coding: utf-8 -*-
"""Tests of the calls report of ``variant_filtration`` with ``tool: varlociraptor``."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.workflows.variant_filtration import VarlociraptorCallsReportStepPart


def _part():
    part = object.__new__(VarlociraptorCallsReportStepPart)
    part.parent = SimpleNamespace(task_name="fdr")
    return part


def test_the_report_reads_the_filtered_vcf_of_the_entity():
    part = _part()
    wildcards = SimpleNamespace(tumor_library="P001-N1-DNA1-WES1")

    assert part._get_input_files_run(wildcards) == {
        "vcf": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz",
        "vcf_tbi": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz.tbi",
    }
    assert part._get_params_run(wildcards) == {"title": "fdr: calls of P001-N1-DNA1-WES1"}


def test_the_report_writes_the_calls_table_and_the_site():
    outputs = _part()._get_output_files_run()

    assert outputs["tsv"] == "work/{tumor_library}/report/calls/{tumor_library}.calls.tsv"
    assert outputs["tsv_md5"] == outputs["tsv"] + ".md5"
    assert is_flagged(outputs["html"], "directory")
    assert str(outputs["html"]) == "work/{tumor_library}/report/calls/html"
    assert VarlociraptorCallsReportStepPart.result_keys == ("tsv", "tsv_md5")
    spec = VarlociraptorCallsReportStepPart.report_outputs["run"]["html"]
    assert (spec.subcategory, spec.htmlindex) == ("varlociraptor calls", "index.html")
