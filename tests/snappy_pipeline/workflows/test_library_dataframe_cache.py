# -*- coding: utf-8 -*-
"""Tests that ``BaseStep`` builds the library table once per run and its dataframe once per
instance."""

from types import SimpleNamespace

import pandas as pd
import pytest

from snappy_pipeline.orchestration import resolve_group_by
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
            "role": ["normal", "normal", "normal"],
            "is_primary": [True, True, True],
            "is_index": [True, False, True],
            "family_id": ["family1", "family1", "family2"],
            "donor_name": ["index", "father", "single"],
            "sex": ["female", "male", "male"],
        }
    )


def _make_step():
    # Bypass ``BaseStep.__init__``, which needs a full project; set only what the
    # dataframe-based properties read.
    step = object.__new__(_DummyStep)
    step.project = SimpleNamespace(library_tables={}, group_by={"dummy": "related"}, sheets=None)
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
        assert step.pedigrees() == {
            "family1": ["index-DNA1", "father-DNA1"],
            "family2": ["single-DNA1"],
        }
        assert step.member_libraries("family1") == ["index-DNA1", "father-DNA1"]
        assert step.primary_library("family1") == "index-DNA1"
        with pytest.raises(ValueError, match="'father-DNA1' is no entity of task 'dummy'"):
            step.member_libraries("father-DNA1")

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
        relationships={"index_lib": {"via": "family_id", "target": "is_index"}},
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
            "role": ["normal", "normal"],
            "is_primary": [True, True],
            "is_index": [True, False],
            "family_id": ["FAM_index", "FAM_index"],
            "donor_name": ["index", "father"],
            "sex": ["female", "male"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    assert step.pedigrees() == {"FAM_index": ["index-DNA1", "father-DNA1"]}
    assert step.pedigree_of("father-DNA1") == "FAM_index"
    assert step.pedigree_index("FAM_index") == "index-DNA1"

    step = _make_step()
    step.config.library_selection = "is_index"
    assert step.pedigrees() == {"FAM_index": ["index-DNA1"]}


def test_a_pedigree_whose_index_has_no_dna_library_has_no_index(mocker):
    # The index has an RNA library only, so it is a group of its own
    df = pd.DataFrame(
        {
            "library_name": ["index-RNA1", "father-DNA1", "mother-DNA1"],
            "extraction_type": ["rna", "dna", "dna"],
            "role": ["normal", "normal", "normal"],
            "is_primary": [True, True, True],
            "is_index": [True, False, False],
            "family_id": ["FAM_index", "FAM_index", "FAM_index"],
            "donor_name": ["index", "father", "mother"],
            "sex": ["female", "male", "female"],
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
            "role": ["normal", "normal", "normal", "normal", "tumor"],
            "is_primary": [True, True, False, True, True],
            "is_index": [False, True, True, True, True],
            "family_id": ["FAM_index"] * 4 + [None],
            "donor_name": ["father", "index", "index", "index", "P1"],
            "sex": ["male", "female", "female", "female", "unknown"],
        }
    )
    mocker.patch.object(abstract, "build_library_dataframe", return_value=df)
    step = _make_step()
    step.config.group_by, step.project.group_by = "family", {"dummy": "family"}
    step.config.library_selection = "family_id == 'FAM_index'"
    assert step.pedigrees() == {"FAM_index": ["father-DNA1", "index-DNA1"]}
    assert step.member_libraries("FAM_index") == ["index-DNA1", "father-DNA1"]
    # A cancer donor without a normal stands for its primary tumor
    step = _make_step()
    step.config.group_by, step.project.group_by = "donor", {"dummy": "donor"}
    step.config.library_selection = "donor_name == 'P1'"
    assert step.member_libraries("P1") == ["P1-T1-DNA1"]
    with pytest.raises(ValueError, match="'P1-T1-DNA1' is no entity of task 'dummy'"):
        step.member_libraries("P1-T1-DNA1")


def test_group_by_comes_from_the_task_then_its_step_then_upstream(mocker):
    """``load_project`` resolves the level of every task once; the steps read ``Project.group_by``."""
    mocker.patch.object(abstract, "build_library_dataframe", return_value=_library_dataframe())
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
    tasks = tuple(SimpleNamespace(name=name, step=step) for name, step in steps.items())
    dependencies = {
        "mapping": {},
        "calling": {"alignments": "mapping"},
        "annotation": {"variants": "calling"},
        "single": {"alignments": "mapping"},
        "combine": {"germline_variants": "calling", "somatic_variants": "single"},
        "phasing": {"combined_variants": "combine"},
    }

    group_by = resolve_group_by(tasks, configs, dependencies)

    assert group_by == {
        "mapping": "library",  # nothing declared
        "calling": "related",  # the step's default: joint calling
        "annotation": "related",  # inherited from the calling task
        "single": "library",  # the task's own choice wins
        "combine": "library",  # per tumor, not the germline calls' grouping
        "phasing": "library",  # the tumor tools phase per tumor
    }
    step = _make_step()
    step.project.group_by, step.task_name = group_by, "annotation"
    assert step.effective_group_by == "family"  # related, resolved by the libraries' family ids


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
