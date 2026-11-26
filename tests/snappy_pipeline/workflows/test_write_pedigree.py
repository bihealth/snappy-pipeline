# -*- coding: utf-8 -*-
"""Tests of the PED files: the family, the individuals (DNA libraries or donors), parents, sex
and affected status in the PED codes."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from snappy_pipeline.library_table import build_donor_dataframe, build_library_dataframe
from snappy_pipeline.sample_sheet import read_sheets
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    WritePedigreeSampleNameStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

SHEET = Path(__file__).parents[1] / "fixtures" / "samplesheet_germline.tsv"


class _Step(BaseStep):
    name = "pedigree_writer"
    produces = [DataSignature(DataType.QC)]


@pytest.fixture
def step():
    sheets = read_sheets([str(SHEET)])
    table = build_library_dataframe(sheets)
    # Bypass ``BaseStep.__init__``; the shared library table is already built
    step = object.__new__(_Step)
    step.task_name = "writer"
    step.config = SimpleNamespace(library_selection=None, relationships=None)
    step.project = SimpleNamespace(
        library_tables={False: table, True: table},
        donor_tables={False: build_donor_dataframe(sheets)},
        sheets=sheets,
    )
    step._library_dataframes = {}
    step._primary_dna_libraries = {}
    return step


def _write(part_class, step, tmp_path) -> list[str]:
    part = object.__new__(part_class)
    part.parent = step
    family = step.primary_dna_libraries()["pedigree"].iloc[0]
    output = tmp_path / "pedigree.ped"
    part.write(part._get_params_run(SimpleNamespace(related_group=family)), str(output))
    return output.read_text().splitlines()


def test_individuals_are_the_dna_libraries(step, tmp_path):
    assert _write(WritePedigreeStepPart, step, tmp_path) == [
        "case001subregion\tcase001subregion-N1-DNA1-WES1\tcase001subregionFather-N1-DNA1-WES1\t"
        "case001subregionMother-N1-DNA1-WES1\t1\t2",
        "case001subregion\tcase001subregionFather-N1-DNA1-WES1\t0\t0\t1\t1",
        "case001subregion\tcase001subregionMother-N1-DNA1-WES1\t0\t0\t2\t1",
    ]


def test_external_exports_name_the_donors(step, tmp_path):
    assert _write(WritePedigreeSampleNameStepPart, step, tmp_path) == [
        "case001subregion\tcase001subregion\tcase001subregionFather\tcase001subregionMother\t1\t2",
        "case001subregion\tcase001subregionFather\t0\t0\t1\t1",
        "case001subregion\tcase001subregionMother\t0\t0\t2\t1",
    ]
