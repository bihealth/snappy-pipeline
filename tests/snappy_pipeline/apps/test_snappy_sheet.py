# -*- coding: utf-8 -*-
"""Tests for ``snappy sheet check`` and ``snappy sheet convert``."""

import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main
from snappy_pipeline.apps.snappy_sheet import convert
from snappy_pipeline.library_table import DataSetInfo, build_library_dataframe
from snappy_pipeline.sample_sheet import read_sheets

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tests/snappy_pipeline/fixtures"
PIPELINES = ROOT / ".tests/test-workflow/pipelines"

#: Every biomedsheets sheet of the repository with its data set type and pedigree field
SHEETS = [
    (FIXTURES / "samplesheet.tsv", "matched_cancer", None),
    (FIXTURES / "samplesheet_rna.tsv", "matched_cancer", None),
    (FIXTURES / "samplesheet_somatic_trio.tsv", "matched_cancer", None),
    (FIXTURES / "samplesheet_germline.tsv", "germline_variants", None),
    (PIPELINES / "snappy-cancer_wes/samplesheet.tsv", "matched_cancer", None),
    (PIPELINES / "snappy-umi/samplesheet.tsv", "matched_cancer", None),
    (PIPELINES / "snappy-germline_wes/samplesheet.tsv", "germline_variants", None),
    (PIPELINES / "snappy-germline_pedigree/samplesheet.tsv", "germline_variants", "familyId"),
]


@pytest.mark.filterwarnings("ignore::biomedsheets.shortcuts.base.MissingDataWarning")
@pytest.mark.parametrize("path, sheet_type, pedigree_field", SHEETS, ids=lambda v: str(v)[-40:])
@pytest.mark.parametrize("naming_scheme", ["only_secondary_id", "secondary_id_pk"])
def test_converted_sheets_load_with_the_names_of_the_old_loader(
    tmp_path, path, sheet_type, pedigree_field, naming_scheme
):
    text = convert(str(path), sheet_type, naming_scheme, pedigree_field)
    (tmp_path / "new.tsv").write_text(text)
    sheets = read_sheets([str(tmp_path / "new.tsv")])

    info = DataSetInfo("old", str(path), sheet_type, False, naming_scheme, pedigree_field)
    old = build_library_dataframe([info], include_background=True).set_index("library_name")
    new = {row.values.library_name: row.values for row in sheets.libraries}
    assert set(new) == set(old.index)
    for name, values in new.items():
        assert values.donor_name == old.loc[name, "donor_name"]
        assert values.sample_name == old.loc[name, "sample_name"]
        assert values.library_type == old.loc[name, "library_type"]
        assert values.extraction_type == old.loc[name, "extraction_type"]
        if sheet_type == "matched_cancer":
            assert values.role == old.loc[name, "role"]
        else:
            assert values.family_id == old.loc[name, "related_group"]
            assert (values.is_index is True) == (old.loc[name, "role"] == "index")
    # Every donor of the old sheet, sequenced or not, has rows
    donors = {row.values.donor_name for row in sheets.rows}
    assert donors >= set(old["donor_name"])


def test_converted_cancer_sheet(tmp_path):
    text = convert(
        str(FIXTURES / "samplesheet_rna.tsv"), "matched_cancer", "only_secondary_id", None
    )
    assert text == textwrap.dedent(
        """\
        donor_name\tsample_name\tis_tumor\tlibrary_name\tlibrary_type\tlibrary_kit\tfolder_name
        case001subregion\tcase001subregion-N1\tfalse\tcase001subregion-N1-DNA1-WES1\tWES\tMedExome\tcase001subregion-N1-DNA1-WES1
        case001subregion\tcase001subregion-T1\ttrue\tcase001subregion-T1-DNA1-WES1\tWES\tMedExome\tcase001subregion-T1-DNA1-WES1
        case001subregion\tcase001subregion-T1\ttrue\tcase001subregion-T1-RNA1-mRNA_seq1\tmRNA_seq\t\tcase001subregion-T1-RNA1-mRNA_seq1
        """
    )


def test_converted_pedigree_sheet_marks_the_index_biomedsheets_chose(tmp_path):
    text = convert(
        str(PIPELINES / "snappy-germline_pedigree/samplesheet.tsv"),
        "germline_variants",
        "only_secondary_id",
        "familyId",
    )
    lines = [line.split("\t") for line in text.splitlines()]
    header = lines[0]
    rows = {row[header.index("donor_name")]: dict(zip(header, row)) for row in lines[1:]}
    assert header[:7] == [
        "donor_name",
        "family_id",
        "father",
        "mother",
        "sex",
        "is_affected",
        "is_index",
    ]
    assert rows["fa"]["is_index"] == "true" and rows["idx"]["is_index"] == ""  # biomedsheets' rule
    assert rows["idx"]["father"] == "fa" and rows["idx"]["mother"] == "mo"
    assert rows["mo"]["is_affected"] == "" and rows["mo"]["mother"] == "mgm"
    assert rows["pgm"]["library_name"] == "" and rows["pgm"]["sex"] == "female"
    assert all(row["family_id"] == "PED1" for row in rows.values())


def test_converted_somatic_trio_marks_the_primary_tumor_and_keeps_the_extra_column(tmp_path):
    text = convert(
        str(FIXTURES / "samplesheet_somatic_trio.tsv"), "matched_cancer", "only_secondary_id", None
    )
    lines = [line.split("\t") for line in text.splitlines()]
    header = lines[0]
    rows = {row[header.index("library_name")]: dict(zip(header, row)) for row in lines[1:]}
    assert rows["case001subregion-T1-DNA1-WES1"]["is_primary"] == "true"
    assert rows["case001subregion-T2-DNA1-WES1"]["is_primary"] == "false"
    assert rows["case001subregion-N1-DNA1-WES1"]["is_primary"] == ""
    assert rows["case001subregion-T2-DNA1-WES1"]["timepoint"] == "relapse"


def test_sheet_check_and_convert_commands(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "sheet",
            "convert",
            str(FIXTURES / "samplesheet.tsv"),
            "--type",
            "matched_cancer",
            "-o",
            str(tmp_path / "s.tsv"),
        ],
    )
    assert result.exit_code == 0, result.output
    result = runner.invoke(main, ["sheet", "check", str(tmp_path / "s.tsv")])
    assert result.exit_code == 0, result.output
    assert f"2 libraries of 1 donors in {tmp_path}/s.tsv" in result.output

    (tmp_path / "bad.tsv").write_text("donor_name\tlibrary_name\tlibrary_type\nP1\tL1\tRNA\n")
    result = runner.invoke(main, ["sheet", "check", str(tmp_path / "bad.tsv")])
    assert result.exit_code == 1
    assert "bad.tsv:2: library_type: 'RNA' is no known value" in result.output


def test_sheet_check_lists_synonyms_and_fixes_them(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("donor_name\tsex\tlibrary_name\tlibrary_type\nP1\tF\tL1\texome\n")
    result = CliRunner().invoke(main, ["sheet", "check", str(path)])
    assert result.exit_code == 0, result.output
    assert f"{path}:2: sex: 'F' is read as 'female'" in result.output
    assert (
        "2 cells are synonyms of their value; snappy sheet check --fix rewrites them"
        in result.output
    )
    result = CliRunner().invoke(main, ["sheet", "check", "--fix", str(path)])
    assert result.exit_code == 0, result.output
    assert f"Rewrote 2 cells of {path}" in result.output
    assert path.read_text() == "donor_name\tsex\tlibrary_name\tlibrary_type\nP1\tfemale\tL1\tWES\n"
