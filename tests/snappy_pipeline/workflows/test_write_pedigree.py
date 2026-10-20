# -*- coding: utf-8 -*-
"""Tests that the PED files are written as biomedsheets wrote them on main."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from biomedsheets.io_tsv import read_germline_tsv_sheet
from biomedsheets.shortcuts import GermlineCaseSheet
from biomedsheets.shortcuts.germline import write_pedigree_to_ped

from snappy_pipeline.workflows.abstract import (
    BaseStep,
    WritePedigreeSampleNameStepPart,
    WritePedigreeStepPart,
    build_library_dataframe,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

SHEET = Path(__file__).parents[1] / "fixtures" / "samplesheet_germline.tsv"


class _Step(BaseStep):
    name = "pedigree_writer"
    produces = [DataSignature(DataType.QC)]


@pytest.fixture
def step():
    with SHEET.open() as f:
        sheet = read_germline_tsv_sheet(f)
    info = SimpleNamespace(
        name="first",
        sheet_type="germline_variants",
        is_background=False,
        pedigree_field_kwargs=None,
    )
    table = build_library_dataframe([info], [sheet])
    # Bypass ``BaseStep.__init__``; the shared library table is already built
    step = object.__new__(_Step)
    step.task_name = "writer"
    step.config = SimpleNamespace(library_selection=None, relationships=None)
    step.project = SimpleNamespace(library_tables={False: table, True: table})
    step._library_dataframes = {}
    step._primary_dna_libraries = {}
    step.sheet = sheet
    return step


def _write(part_class, step, tmp_path) -> list[str]:
    part = object.__new__(part_class)
    part.parent = step
    index = step.primary_dna_libraries()["pedigree"].iloc[0]
    output = tmp_path / "pedigree.ped"
    part.run(SimpleNamespace(index_ngs_library=index), str(output))
    return output.read_text().splitlines()


def test_individuals_are_the_dna_libraries_as_on_main(step, tmp_path):
    expected = tmp_path / "main.ped"
    write_pedigree_to_ped(GermlineCaseSheet(step.sheet).cohort.pedigrees[0], str(expected))
    assert _write(WritePedigreeStepPart, step, tmp_path) == expected.read_text().splitlines()


def test_external_exports_name_the_donors(step, tmp_path):
    family = "FAM_case001subregion-000001"
    assert _write(WritePedigreeSampleNameStepPart, step, tmp_path) == [
        f"{family}\tcase001subregion-000001\tcase001subregionFather-000005\t"
        "case001subregionMother-000009\t1\t2",
        f"{family}\tcase001subregionFather-000005\t0\t0\t1\t1",
        f"{family}\tcase001subregionMother-000009\t0\t0\t2\t1",
    ]
