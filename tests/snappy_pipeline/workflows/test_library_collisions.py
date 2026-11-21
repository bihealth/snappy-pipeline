# -*- coding: utf-8 -*-
"""Tests for libraries of the same name in several data sets."""

import pandas as pd
import pytest

from snappy_pipeline import orchestration
from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.library_table import select_libraries

FIX = "library_selection: \"data_set == 'batch1'\""


def _rows(*data_sets, name="P001-N1-DNA1-WES1"):
    return [
        {
            "library_name": name,
            "data_set": data_set,
            "folder_name": name,
            "extraction_type": "dna",
            "kind": "germline",
        }
        for data_set in data_sets
    ]


def _step(task_name, rows, selected):
    # library_row only reads these attributes, so skip the full step setup
    step = object.__new__(BaseStep)
    step.task_name = task_name
    step.__dict__["_library_rows"] = {row["library_name"]: [] for row in rows}
    for row in rows:
        step.__dict__["_library_rows"][row["library_name"]].append(row)
    step.__dict__["_selected_libraries"] = {(r["data_set"], r["library_name"]) for r in selected}
    return step


def test_a_selection_with_colliding_names_fails_with_the_fix():
    df = pd.DataFrame(_rows("batch1", "batch2"))
    with pytest.raises(
        ConfigError, match="several data sets: P001-N1-DNA1-WES1 \\(batch1, batch2\\)"
    ) as error:
        select_libraries(df, None)
    assert FIX in str(error.value)
    assert len(select_libraries(df, "data_set == 'batch2'")) == 1


def test_a_related_group_in_several_data_sets_fails_with_the_fix():
    rows = _rows("batch1") + _rows("batch2", name="P001-N2-DNA1-WES1")
    for row in rows:
        row["related_group"] = "FAM_P001"
    with pytest.raises(
        ConfigError, match="related groups of the same name in several data sets: FAM_P001"
    ) as error:
        select_libraries(pd.DataFrame(rows), None)
    assert FIX in str(error.value)
    assert len(select_libraries(pd.DataFrame(rows), "data_set == 'batch2'")) == 1


def test_library_selection_decides_between_libraries_of_one_name():
    rows = _rows("batch1", "batch2") + _rows("batch1", name="P002-N1-DNA1-WES1")
    step = _step("mapping", rows, selected=rows[1:])
    assert step.library_row("P001-N1-DNA1-WES1")["data_set"] == "batch2"
    # An unambiguous name resolves whatever the selection keeps
    assert step.library_row("P002-N1-DNA1-WES1")["data_set"] == "batch1"

    step = _step("mapping", rows, selected=rows)
    with pytest.raises(
        ConfigError, match="^libraries of the same name in several data sets"
    ) as error:
        step.library_row("P001-N1-DNA1-WES1")
    assert FIX in str(error.value)


def test_upstream_task_must_write_a_result_for_a_related_group(monkeypatch):
    rows = _rows("batch1")
    mapping = _step("mapping", rows, selected=rows)
    mapping.__dict__["output_entities"] = ["P001-N1-DNA1-WES1"]
    mapping.__dict__["effective_group_by"] = "library"
    calling = _step("calling", rows, selected=rows)
    monkeypatch.setattr(orchestration, "_TASK_INSTANCES", {"mapping": mapping, "calling": calling})
    with pytest.raises(
        ConfigError,
        match="task 'mapping' writes one result per library, but 'FAM_P001' is a related group; "
        "set group_by: library on task 'calling', or group_by: related on 'mapping'",
    ):
        calling._check_upstream_library("mapping", "FAM_P001")

    mapping.__dict__["effective_group_by"] = "related"
    mapping.__dict__["output_entities"] = ["FAM_P002"]
    with pytest.raises(ConfigError, match="writes no result for the related group 'FAM_P001'"):
        calling._check_upstream_library("mapping", "FAM_P001")

    mapping.__dict__["output_entities"] = ["FAM_P001"]
    calling._check_upstream_library("mapping", "FAM_P001")


def test_upstream_task_must_have_the_library_of_the_same_data_set(monkeypatch):
    rows = _rows("batch1", "batch2")
    mapping = _step("mapping", rows, selected=rows[:1])
    calling = _step("calling", rows, selected=rows[1:])
    monkeypatch.setattr(orchestration, "_TASK_INSTANCES", {"mapping": mapping, "calling": calling})
    with pytest.raises(
        ConfigError, match="from data set 'batch2', but the upstream task 'mapping'"
    ):
        calling._check_upstream_library("mapping", "P001-N1-DNA1-WES1")

    calling.__dict__["_selected_libraries"] = mapping._selected_libraries
    calling._check_upstream_library("mapping", "P001-N1-DNA1-WES1")
