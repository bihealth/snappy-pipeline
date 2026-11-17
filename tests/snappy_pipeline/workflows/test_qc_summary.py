# -*- coding: utf-8 -*-
"""Tests of the QC summary step part: the files it reads per entity, its labels, outputs and
result files, and the subclasses of the steps."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.ngs_data_qc import PicardSummaryStepPart
from snappy_pipeline.workflows.ngs_mapping import BamQcSummaryStepPart
from snappy_pipeline.workflows.variant_calling import BcftoolsStatsSummaryStepPart


class _Source:
    """A step part that writes two files per entity, with the result wildcards it is given."""

    def __init__(self, wildcards):
        self.wildcards = wildcards

    def get_output_files(self, action):
        assert action == "run"
        return {
            "txt": "work/{index,[^/]+}/report/{index}.{donor}.txt",
            "txt_md5": "work/{index}/report/{index}.{donor}.txt.md5",
        }

    def _result_wildcards(self, entity):
        return {"index": [entity], **self.wildcards.get(entity, {})}


class _Summary(QcSummaryStepPart):
    source = ("stats", "run")
    source_keys = {"stats": "txt"}
    format = "stats"
    report_outputs = {"run": {"html": SimpleNamespace(subcategory="Statistics")}}


def _part(entities, wildcards=None, label_wildcard=None):
    part = object.__new__(_Summary)
    part.parent = SimpleNamespace(
        output_entities=entities, sub_steps={"stats": _Source(wildcards or {})}
    )
    part.label_wildcard = label_wildcard
    return part


def test_inputs_are_the_source_files_of_every_entity_in_order():
    part = _part(["b", "a"], {"b": {"donor": ["b"]}, "a": {"donor": ["a"]}})

    inputs = part._get_input_files_run(None)
    params = part._get_params_run(None)

    assert inputs == {"stats": ["work/b/report/b.b.txt", "work/a/report/a.a.txt"]}
    assert params == {"format": "stats", "labels": ["b", "a"], "title": "Statistics"}


def test_a_label_wildcard_labels_every_file_of_an_entity():
    part = _part(["idx"], {"idx": {"donor": ["idx", "father"]}}, label_wildcard="donor")

    assert part._get_input_files_run(None) == {
        "stats": ["work/idx/report/idx.idx.txt", "work/idx/report/idx.father.txt"]
    }
    assert part._get_params_run(None)["labels"] == ["idx", "father"]


def test_outputs_are_the_summary_tsv_and_the_site_directory():
    part = _part(["a"], {"a": {"donor": ["a"]}})

    outputs = part._get_output_files_run()

    assert outputs["tsv"] == "work/qc_summary/out/qc_summary.tsv"
    assert outputs["tsv_md5"] == "work/qc_summary/out/qc_summary.tsv.md5"
    assert is_flagged(outputs["html"], "directory")
    assert str(outputs["html"]) == "work/qc_summary/report/html"
    assert part._get_log_file_run() == {"log": "work/qc_summary/log/qc_summary.log"}


def test_result_files_are_the_tsv_once_or_nothing_without_entities():
    assert _part(["a", "b"], {"a": {"donor": ["a"]}, "b": {"donor": ["b"]}}).get_result_files() == [
        "output/qc_summary/out/qc_summary.tsv",
        "output/qc_summary/out/qc_summary.tsv.md5",
    ]
    assert _part([]).get_result_files() == []


def test_the_steps_name_their_sources():
    mapping = object.__new__(BamQcSummaryStepPart)
    mapping.config = SimpleNamespace(tool="bwa_mem2")
    assert mapping.source == ("bwa_mem2", "run")
    assert mapping.format == "samtools"
    assert set(mapping.source_keys) == {"bamstats", "flagstats", "idxstats"}

    picard = object.__new__(PicardSummaryStepPart)
    picard.config = SimpleNamespace(
        picard=SimpleNamespace(programs=["CollectHsMetrics", "CollectInsertSizeMetrics"])
    )
    assert picard.source == ("picard", "metrics")
    assert picard.source_keys == {
        "CollectHsMetrics": "CollectHsMetrics",
        "CollectInsertSizeMetrics": "CollectInsertSizeMetrics",
    }

    assert BcftoolsStatsSummaryStepPart.label_wildcard == "donor_library_name"
    assert BcftoolsStatsSummaryStepPart.source == ("bcftools_stats", "run")
