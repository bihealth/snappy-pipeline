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


def test_the_fingerprint_checks_pass_the_sheet_facts_as_options():
    import pandas as pd

    from snappy_pipeline.workflows.ngs_data_qc import NgsChewSummaryStepPart

    part = object.__new__(NgsChewSummaryStepPart)
    table = pd.DataFrame(
        [
            {
                "library_name": "idx-N1-DNA1-WES1",
                "donor_name": "idx",
                "related_group": "FAM_idx",
                "sex": "male",
            },
            {
                "library_name": "mo-N1-DNA1-WES1",
                "donor_name": "mo",
                "related_group": "FAM_idx",
                "sex": "female",
            },
        ]
    )
    donors = pd.DataFrame(
        [
            {
                "donor_name": "fa",
                "father_name": "0",
                "mother_name": "0",
                "sex": "male",
                "library_name": "",
            },
            {
                "donor_name": "mo",
                "father_name": "0",
                "mother_name": "0",
                "sex": "female",
                "library_name": "mo-N1-DNA1-WES1",
            },
            {
                "donor_name": "idx",
                "father_name": "fa",
                "mother_name": "mo",
                "sex": "male",
                "library_name": "idx-N1-DNA1-WES1",
            },
        ]
    )
    part.parent = SimpleNamespace(
        build_library_dataframe=lambda: table,
        dna_libraries=lambda: ["idx-N1-DNA1-WES1", "mo-N1-DNA1-WES1"],
        pedigree_donors=lambda pedigree: donors,
        output_entities=["idx-N1-DNA1-WES1", "mo-N1-DNA1-WES1"],
        sub_steps={},
    )

    assert NgsChewSummaryStepPart.source == ("ngs_chew", "compare")
    assert set(NgsChewSummaryStepPart.source_keys) == {"relatedness", "ibs0", "mask", "stats"}
    assert part.entities() == ["cohort"]
    options = part.summary_options()
    assert options["libraries"] == [
        {"library": "idx-N1-DNA1-WES1", "donor": "idx", "pedigree": "FAM_idx", "sex": "male"},
        {"library": "mo-N1-DNA1-WES1", "donor": "mo", "pedigree": "FAM_idx", "sex": "female"},
    ]
    assert options["donors"][0] == {
        "donor": "fa",
        "pedigree": "FAM_idx",
        "father": "0",
        "mother": "0",
        "sex": "male",
        "library": "",
    }
    assert len(options["donors"]) == 3


def test_summary_options_reach_the_params():
    class _WithOptions(_Summary):
        def summary_options(self):
            return {"libraries": [{"library": "a"}]}

    part = object.__new__(_WithOptions)
    part.parent = SimpleNamespace(
        output_entities=["a"], sub_steps={"stats": _Source({"a": {"donor": ["a"]}})}
    )
    part.label_wildcard = None

    assert part._get_params_run(None)["libraries"] == [{"library": "a"}]
    assert _part(["a"], {"a": {"donor": ["a"]}})._get_params_run(None) == {
        "format": "stats",
        "labels": ["a"],
        "title": "Statistics",
    }
