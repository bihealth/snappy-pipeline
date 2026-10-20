# -*- coding: utf-8 -*-
"""Tests that ``BaseStep`` builds the library table once per run and its dataframe once per
instance."""

from types import SimpleNamespace

import pandas as pd

from snappy_pipeline.workflows import abstract
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType


class _DummyStep(BaseStep):
    name = "dummy"
    produces = [DataSignature(DataType.QC)]


def _library_dataframe():
    return pd.DataFrame(
        {
            "library_name": ["index-DNA1", "father-DNA1", "single-DNA1"],
            "extraction_type": ["dna", "dna", "dna"],
            "kind": ["germline", "germline", "germline"],
            "role": ["index", "father", "index"],
            "is_primary": [True, False, True],
            "cohort_name": ["family1", "family1", "family2"],
        }
    )


def _make_step():
    # Bypass ``BaseStep.__init__``, which needs a full project; set only what the
    # dataframe-based properties read.
    step = object.__new__(_DummyStep)
    step.data_set_infos = []
    step.sheets = []
    step.project = SimpleNamespace(library_tables={})
    step.config = SimpleNamespace(library_selection=None, group_by="cohort", relationships=None)
    step.task_name = "dummy"
    step._library_dataframes = {}
    step._primary_dna_libraries = {}
    return step


def test_library_dataframe_is_built_once(mocker):
    build = mocker.patch.object(
        abstract, "build_library_dataframe", return_value=_library_dataframe()
    )
    step = _make_step()

    for _ in range(3):
        step.build_library_dataframe()
        assert step.output_entities == ["index-DNA1", "single-DNA1"]
        assert step.cohort_members == {
            "family1": ["index-DNA1", "father-DNA1"],
            "family2": ["single-DNA1"],
        }
        assert step.get_cohort_libraries("father-DNA1") == ["index-DNA1", "father-DNA1"]
        assert step.get_cohort_libraries("unknown-DNA1") == ["unknown-DNA1"]

    assert build.call_count == 1


def test_tasks_share_the_library_table_of_the_run(mocker):
    build = mocker.patch.object(
        abstract, "build_library_dataframe", return_value=_library_dataframe()
    )
    first, second = _make_step(), _make_step()
    second.project = first.project
    second.config = SimpleNamespace(
        library_selection=None,
        group_by=None,
        relationships={"index_lib": {"via": "cohort_name", "target": "role == 'index'"}},
    )

    assert "index_lib" not in first.build_library_dataframe().columns
    assert second.build_library_dataframe()["index_lib"].tolist() == [
        "index-DNA1",
        "index-DNA1",
        "single-DNA1",
    ]
    assert build.call_count == 1


def test_library_selection_narrows_the_pedigrees(mocker):
    df = pd.DataFrame(
        {
            "library_name": ["index-DNA1", "father-DNA1"],
            "extraction_type": ["dna", "dna"],
            "kind": ["germline", "germline"],
            "role": ["index", "father"],
            "cohort_name": ["index-DNA1", "index-DNA1"],
            "donor_name": ["index", "father"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    assert step.pedigrees() == {"index-DNA1": ["index-DNA1", "father-DNA1"]}

    step = _make_step()
    step.config.library_selection = "role != 'father'"
    assert step.pedigrees() == {"index-DNA1": ["index-DNA1"]}


def test_group_by_comes_from_the_task_then_its_step_then_upstream():
    configs = {
        "mapping": SimpleNamespace(group_by=None, tool="bwa"),
        "calling": SimpleNamespace(group_by=None, tool="gatk4_hc_joint"),
        "annotation": SimpleNamespace(group_by=None, tool="vep"),
        "single": SimpleNamespace(group_by="library", tool="gatk4_hc_joint"),
    }
    steps = {"mapping": "ngs_mapping", "calling": "variant_calling"}
    steps |= {"annotation": "variant_annotation", "single": "variant_calling"}
    project = SimpleNamespace(
        task_configs=configs,
        dependencies={
            "calling": {"alignments": "mapping"},
            "annotation": {"variants": "calling"},
            "single": {"alignments": "mapping"},
        },
        task=lambda name: SimpleNamespace(step=steps[name]),
    )

    def effective(name):
        step = object.__new__(_DummyStep)
        step.config, step.project, step.task_name = configs[name], project, name
        return step.effective_group_by

    assert effective("mapping") == "library"  # nothing declared
    assert effective("calling") == "cohort"  # the step's default: joint calling
    assert effective("annotation") == "cohort"  # inherited from the calling task
    assert effective("single") == "library"  # the task's own choice wins


def test_background_selection_narrows_only_the_background_libraries():
    table = pd.DataFrame(
        {
            "library_name": ["case-DNA1", "bg1-DNA1", "bg2-DNA1"],
            "data_set": ["cases", "cohort_a", "cohort_b"],
            "is_background": [False, True, True],
        }
    )
    step = _make_step()
    step.project.library_tables[True] = table
    step.config.background_selection = "data_set == 'cohort_b'"
    assert step.build_library_dataframe(include_background=True)["library_name"].tolist() == [
        "case-DNA1",
        "bg2-DNA1",
    ]
