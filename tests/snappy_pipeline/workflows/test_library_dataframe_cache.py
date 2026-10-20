# -*- coding: utf-8 -*-
"""Tests that ``BaseStep`` builds the library dataframe only once per instance."""

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
    step.shortcut_sheets = []
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
