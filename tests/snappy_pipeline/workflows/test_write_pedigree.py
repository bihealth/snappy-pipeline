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


#: Two cousins (c1 the index, c2) through the siblings pa and pb, whose parents have no data
COUSINS = """\
family_id\tdonor_name\tfather\tmother\tsex\tis_affected\tis_index\tsample_name\tlibrary_name\tlibrary_type
FAM\tgf\t\t\tmale\t\t\t\t\t
FAM\tgm\t\t\tfemale\t\t\t\t\t
FAM\tpa\tgf\tgm\tmale\tfalse\t\tpa-N1\tpa-N1-DNA1-WES1\tWES
FAM\tpb\tgf\tgm\tfemale\tfalse\t\tpb-N1\tpb-N1-DNA1-WES1\tWES
FAM\tc1\tpa\t\tmale\ttrue\ttrue\tc1-N1\tc1-N1-DNA1-WES1\tWES
FAM\tc2\t\tpb\tfemale\ttrue\t\tc2-N1\tc2-N1-DNA1-WES1\tWES
"""


@pytest.fixture
def step():
    return _make_step(read_sheets([str(SHEET)]))


def _make_step(sheets):
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


def test_relatives_through_donors_without_data_are_unrelated_in_the_ped(tmp_path):
    """The grandparents connect the cousins with donor-only rows; the PED names the donors with a
    DNA library only, so the siblings pa and pb have no parents there."""
    (tmp_path / "cousins.tsv").write_text(COUSINS)
    sheets = read_sheets([str(tmp_path / "cousins.tsv")])
    step = _make_step(sheets)
    donors = step.pedigree_donors("FAM")
    assert donors["donor_name"].tolist() == ["gf", "gm", "pa", "pb", "c1", "c2"]
    assert donors["father_name"].tolist() == ["0", "0", "gf", "gf", "pa", "0"]
    assert step.pedigree_index("FAM") == "c1-N1-DNA1-WES1"
    assert _write(WritePedigreeStepPart, step, tmp_path) == [
        "FAM\tpa-N1-DNA1-WES1\t0\t0\t1\t1",
        "FAM\tpb-N1-DNA1-WES1\t0\t0\t2\t1",
        "FAM\tc1-N1-DNA1-WES1\tpa-N1-DNA1-WES1\t0\t1\t2",
        "FAM\tc2-N1-DNA1-WES1\t0\tpb-N1-DNA1-WES1\t2\t2",
    ]
