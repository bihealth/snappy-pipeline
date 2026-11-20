# -*- coding: utf-8 -*-
"""Tests that ``BaseStep`` builds the library table once per run and its dataframe once per
instance."""

from types import SimpleNamespace

import pandas as pd
import pytest

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
            "related_group": ["family1", "family1", "family2"],
        }
    )


def _make_step():
    # Bypass ``BaseStep.__init__``, which needs a full project; set only what the
    # dataframe-based properties read.
    step = object.__new__(_DummyStep)
    step.data_set_infos = []
    step.sheets = []
    step.project = SimpleNamespace(library_tables={})
    step.config = SimpleNamespace(library_selection=None, group_by="related", relationships=None)
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
        assert step.output_entities == ["family1", "family2"]
        assert step.related_libraries == {
            "family1": ["index-DNA1", "father-DNA1"],
            "family2": ["single-DNA1"],
        }
        assert step.get_related_libraries("family1") == ["index-DNA1", "father-DNA1"]
        assert step.primary_library("family1") == "index-DNA1"
        with pytest.raises(ValueError, match="'father-DNA1' is no entity of task 'dummy'"):
            step.get_related_libraries("father-DNA1")

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
        relationships={"index_lib": {"via": "related_group", "target": "role == 'index'"}},
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
            "is_primary": [True, False],
            "related_group": ["FAM_index", "FAM_index"],
            "donor_name": ["index", "father"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    assert step.pedigrees() == {"FAM_index": ["index-DNA1", "father-DNA1"]}
    assert step.pedigree_of("father-DNA1") == "FAM_index"
    assert step.pedigree_index("FAM_index") == "index-DNA1"

    step = _make_step()
    step.config.library_selection = "role != 'father'"
    assert step.pedigrees() == {"FAM_index": ["index-DNA1"]}


def test_a_pedigree_whose_index_has_no_dna_library_has_no_index(mocker):
    # The index has an RNA library only, so it is a group of its own
    df = pd.DataFrame(
        {
            "library_name": ["index-RNA1", "father-DNA1", "mother-DNA1"],
            "extraction_type": ["rna", "dna", "dna"],
            "kind": ["germline", "germline", "germline"],
            "role": ["index", "father", "mother"],
            "is_primary": [False, False, False],
            "related_group": ["FAM_index", "FAM_index", "FAM_index"],
            "donor_name": ["index", "father", "mother"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    assert step.pedigrees() == {"FAM_index": ["father-DNA1", "mother-DNA1"]}
    assert step.pedigree_index("FAM_index") is None
    assert step.member_libraries("FAM_index") == ["father-DNA1", "mother-DNA1"]


def test_a_joint_call_takes_one_dna_library_per_donor_the_index_first(mocker):
    df = pd.DataFrame(
        {
            "library_name": ["father-DNA1", "index-DNA1", "index-DNA2", "index-RNA1", "P1-T1-DNA1"],
            "extraction_type": ["dna", "dna", "dna", "rna", "dna"],
            "kind": ["germline"] * 4 + ["cancer"],
            "role": ["father", "index", "index", "index", "tumor"],
            "is_primary": [False, True, False, False, True],
            "related_group": ["FAM_index"] * 4 + ["PAT_P1"],
            "donor_name": ["father", "index", "index", "index", "P1"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    assert step.get_related_libraries("FAM_index") == df["library_name"].tolist()[:4]
    assert step.pedigrees() == {"FAM_index": ["father-DNA1", "index-DNA1"]}
    assert step.member_libraries("FAM_index") == ["index-DNA1", "father-DNA1"]
    # A cancer donor without a normal stands for its primary tumor
    assert step.member_libraries("PAT_P1") == ["P1-T1-DNA1"]
    with pytest.raises(ValueError, match="'P1-T1-DNA1' is no entity of task 'dummy'"):
        step.member_libraries("P1-T1-DNA1")


def test_group_by_comes_from_the_task_then_its_step_then_upstream():
    configs = {
        "mapping": SimpleNamespace(group_by=None, tool="bwa"),
        "calling": SimpleNamespace(group_by=None, tool="gatk4_hc_joint"),
        "annotation": SimpleNamespace(group_by=None, tool="vep"),
        "single": SimpleNamespace(group_by="library", tool="gatk4_hc_joint"),
        "combine": SimpleNamespace(group_by=None),
        "phasing": SimpleNamespace(group_by=None, tool="whatshap_tumor"),
    }
    steps = {"mapping": "ngs_mapping", "calling": "variant_calling"}
    steps |= {"annotation": "variant_annotation", "single": "variant_calling"}
    steps |= {"combine": "combine_variants", "phasing": "variant_phasing"}
    project = SimpleNamespace(
        task_configs=configs,
        dependencies={
            "calling": {"alignments": "mapping"},
            "annotation": {"variants": "calling"},
            "single": {"alignments": "mapping"},
            "combine": {"germline_variants": "calling", "somatic_variants": "single"},
            "phasing": {"combined_variants": "combine"},
        },
        task=lambda name: SimpleNamespace(step=steps[name]),
    )

    def effective(name):
        step = object.__new__(_DummyStep)
        step.config, step.project, step.task_name = configs[name], project, name
        return step.effective_group_by

    assert effective("mapping") == "library"  # nothing declared
    assert effective("calling") == "related"  # the step's default: joint calling
    assert effective("annotation") == "related"  # inherited from the calling task
    assert effective("single") == "library"  # the task's own choice wins
    assert effective("combine") == "library"  # per tumor, not the germline calls' grouping
    assert effective("phasing") == "library"  # the tumor tools phase per tumor


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
