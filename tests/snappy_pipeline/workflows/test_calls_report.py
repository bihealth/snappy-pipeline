# -*- coding: utf-8 -*-
"""Tests of the shared calls report step part of the VCF steps."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart


def _paths(config, library_name):
    return {
        "vcf": f"work/{library_name}/out/{library_name}.vcf.gz",
        "vcf_tbi": f"work/{library_name}/out/{library_name}.vcf.gz.tbi",
    }


def _part(tags=frozenset({"germline", "filtered"}), **settings):
    part = object.__new__(CallsReportStepPart)
    signatures = {"fdr": (DataSignature(DataType.VARIANTS, frozenset(tags)),)}
    part.parent = SimpleNamespace(
        task_name="fdr",
        config=SimpleNamespace(calls_report=CallsReport(**settings)),
        get_output_paths=_paths,
        project=SimpleNamespace(signatures=signatures),
    )
    return part


def test_the_report_reads_the_tasks_own_vcf_of_the_entity():
    part = _part(page_size=10)
    wildcards = SimpleNamespace(library_name="P001-N1-DNA1-WES1")

    assert part._get_input_files_run(wildcards) == {
        "vcf": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz",
        "vcf_tbi": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz.tbi",
    }
    assert part._get_params_run(wildcards) == {
        "title": "fdr: calls of P001-N1-DNA1-WES1",
        "page_size": 10,
        "max_in_memory_rows": 50000,
        "linkouts": [
            "gnomad",
            "clinvar",
            "dbsnp",
            "ensembl",
            "genecards",
            "omim",
            "varlociraptor_inspect",
        ],
        "assembly": None,
    }


def test_the_assembly_is_the_tag_of_the_tasks_own_variants():
    tagged = _part(tags={"germline", "annotated", "mehari", "GRCh37"}, linkouts=["gnomad"])

    assert tagged.assembly() == "GRCh37"
    params = tagged._get_params_run(SimpleNamespace(library_name="FAM_P1"))
    assert (params["assembly"], params["linkouts"]) == ("GRCh37", ["gnomad"])
    assert _part().assembly() is None


def test_the_report_writes_the_calls_table_and_the_site():
    outputs = _part()._get_output_files_run()

    assert outputs["tsv"] == "work/{library_name}/report/calls/{library_name}.calls.tsv"
    assert outputs["tsv_md5"] == outputs["tsv"] + ".md5"
    assert is_flagged(outputs["html"], "directory")
    assert str(outputs["html"]) == "work/{library_name}/report/calls/html"
    assert CallsReportStepPart.result_keys == ("tsv", "tsv_md5")
    spec = CallsReportStepPart.report_outputs["run"]["html"]
    assert (spec.subcategory, spec.htmlindex) == ("calls", "index.html")
    assert _part()._get_log_file_run() == {
        "log": "work/{library_name}/log/{library_name}.calls_report.log"
    }
