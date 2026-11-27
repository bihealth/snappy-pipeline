# -*- coding: utf-8 -*-
"""Tests for the plain TSV sample sheets (``snappy_pipeline.sample_sheet``)."""

import pytest

from snappy_pipeline.sample_sheet import (
    COLUMNS,
    Level,
    SheetError,
    header,
    primary_library,
    read_sheets,
)

CANCER = """\
donor_name\tsex\tsample_name\tis_tumor\tpurity\tlibrary_name\tlibrary_type\tlibrary_kit\tfolder_name
P001\tfemale\tP001-N1\tfalse\t\tP001-N1-DNA1-WES1\tWES\tMedExome\tP001-N1-DNA1-WES1
P001\tfemale\tP001-T1\ttrue\t0.8\tP001-T1-DNA1-WES1\tWES\tMedExome\tP001-T1-DNA1-WES1
P001\tfemale\tP001-T1\ttrue\t0.8\tP001-T1-RNA1-mRNA_seq1\tmRNA_seq\t\tP001-T1-RNA1-mRNA_seq1
"""

GERMLINE = """\
# a trio with an unsequenced grandmother
family_id\tdonor_name\tfather\tmother\tsex\tis_affected\tis_index\thpo_terms\tsample_name\tlibrary_name\tlibrary_type\tlibrary_kit

FAM1\tidx\tfa\tmo\tmale\ttrue\ttrue\tHP:0000001,HP:0000002\tidx-N1\tidx-N1-DNA1-WES1\tWES\tV6
FAM1\tfa\t\t\tmale\tfalse\t\t\tfa-N1\tfa-N1-DNA1-WES1\tWES\tV6
FAM1\tmo\t\tgm\tfemale\tfalse\t\t\tmo-N1\tmo-N1-DNA1-WES1\tWES\tV6
FAM1\tgm\t\t\tfemale\t\t\t\t\t\t\t
"""


def _write(tmp_path, text, name="samplesheet.tsv"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _problems(tmp_path, text, name="samplesheet.tsv"):
    """The problems of the sheet ``text`` as strings, without the directory of the sheet."""
    with pytest.raises(SheetError) as excinfo:
        read_sheets([_write(tmp_path, text, name)])
    return [str(p).replace(f"{tmp_path}/", "") for p in excinfo.value.problems]


def test_cancer_sheet(tmp_path):
    sheets = read_sheets([_write(tmp_path, CANCER)])
    rows = [row.values for row in sheets.libraries]
    assert [row.library_name for row in rows] == [
        "P001-N1-DNA1-WES1",
        "P001-T1-DNA1-WES1",
        "P001-T1-RNA1-mRNA_seq1",
    ]
    assert [row.role for row in rows] == ["normal", "tumor", "tumor"]
    assert [row.extraction_type for row in rows] == ["dna", "dna", "rna"]
    assert rows[1].purity == 0.8 and rows[0].purity is None and rows[2].library_kit is None
    assert sheets.libraries[1].line == 3
    assert sheets.columns == [
        "donor_name", "sex", "sample_name", "is_tumor", "purity",
        "library_name", "library_type", "library_kit", "folder_name",
    ]  # fmt: skip


def test_germline_sheet_with_comments_blank_lines_and_a_donor_without_libraries(tmp_path):
    sheets = read_sheets([_write(tmp_path, GERMLINE.replace("\n", "\r\n"))])
    assert [row.line for row in sheets.rows] == [4, 5, 6, 7]
    assert len(sheets.libraries) == 3
    gm = sheets.rows[3].values
    assert not gm.is_library and gm.sex == "female" and gm.is_affected is None
    idx = sheets.rows[0].values
    assert idx.is_index is True and idx.hpo_terms == "HP:0000001,HP:0000002" and idx.father == "fa"


def test_mixed_sheet_and_extra_columns(tmp_path):
    text = (
        "family_id\tdonor_name\tsex\tis_index\tis_tumor\tlibrary_name\tlibrary_type\ttimepoint\n"
        "F1\tP1\tmale\ttrue\tfalse\tP1-N1-DNA1-WGS1\tWGS\t\n"
        "F1\tP1\tmale\ttrue\ttrue\tP1-T1-DNA1-WGS1\tWGS\tprimary\n"
        "F1\tP2\tfemale\t\tfalse\tP2-N1-DNA1-WGS1\tWGS\t\n"
    )
    sheets = read_sheets([_write(tmp_path, "﻿" + text)])
    rows = [row.values for row in sheets.rows]
    assert [row.role for row in rows] == ["normal", "tumor", "normal"]
    assert rows[1].extra("timepoint") == "primary" and rows[0].extra("timepoint") is None
    assert sheets.columns[-1] == "timepoint"


def test_several_sheets_are_checked_together(tmp_path):
    a = _write(tmp_path, "donor_name\tlibrary_name\tlibrary_type\nP1\tL1\tWES\n", "a.tsv")
    b = _write(tmp_path, "donor_name\tlibrary_name\tlibrary_type\nP2\tL1\tWGS\n", "b.tsv")
    with pytest.raises(SheetError, match=r"b.tsv:2: library_name: 'L1' is also on line 2"):
        read_sheets([a, b])
    assert len(read_sheets([a]).libraries) == 1


def test_a_missing_sheet(tmp_path):
    with pytest.raises(SheetError, match="none.tsv: cannot read the sample sheet"):
        read_sheets([str(tmp_path / "none.tsv")])


@pytest.mark.parametrize(
    "text, problems",
    [
        ("", ["samplesheet.tsv: the sample sheet has no header line"]),
        ("# only comments\n\n", ["samplesheet.tsv: the sample sheet has no header line"]),
        (
            "donorName\tlibrary_name\tlibrary_type\n",
            [
                "samplesheet.tsv:1: donorName: column names are snake_case: ^[a-z][a-z0-9_]*$",
                "samplesheet.tsv:1: donor_name: the column is required",
            ],
        ),
        (
            "donor_name\tlibrary_name\tlibrary_type\trole\n",
            ["samplesheet.tsv:1: role: the library table derives this column; rename it"],
        ),
        (
            "donor_name\tlibrary_name\tlibrary_type\tlibrary_name\n",
            ["samplesheet.tsv:1: library_name: the column is repeated"],
        ),
        ("donor_name\tlibrary_name\n", ["samplesheet.tsv:1: library_type: the column is required"]),
    ],
)
def test_header_problems(tmp_path, text, problems):
    assert _problems(tmp_path, text) == problems


HEADER = "donor_name\tsex\tsample_name\tis_tumor\tpurity\tlibrary_name\tlibrary_type\tfolder_name\thpo_terms\tread_structure\n"


@pytest.mark.parametrize(
    "row, problems",
    [
        (
            "P1\tmale\tS1\tfalse\t\tL1\tWES\tF1\t",
            [
                "2: 9 cells, the header has 10; a missing value is an empty cell (some editors drop trailing tabs)"
            ],
        ),
        (
            "P1\tmale\tS1\tfalse\t.\tL1\tWES\tF1\t\t",
            ["2: purity: '.' is no value; leave the cell empty"],
        ),
        (
            "P1 \tmale\tS1\tfalse\t\tL1\tWES\tF1\t\t",
            ["2: donor_name: 'P1 ' has leading or trailing whitespace"],
        ),
        ("P1\tmale\tS1\tmaybe\t\tL1\tWES\tF1\t\t", ["2: is_tumor: 'maybe' is no known value"]),
        (
            "P1\tdiverse\tS1\tfalse\t\tL1\tWES\tF1\t\t",
            ["2: sex: 'diverse' is no known value"],
        ),
        (
            "P1\tmale\tS1\ttrue\t1.5\tL1\tWES\tF1\t\t",
            ["2: purity: Input should be less than or equal to 1"],
        ),
        (
            "P1\tmale\tS1\ttrue\t0,8\tL1\tWES\tF1\t\t",
            ["2: purity: Input should be a valid number, unable to parse string as a number"],
        ),
        (
            "P1\tmale\tS1\tfalse\t0.8\tL1\tWES\tF1\t\t",
            ["2: purity: only a tumor sample has a purity"],
        ),
        (
            "P1\tmale\tS1\tfalse\t\tL1\tWES\ta/b\t\t",
            ["2: folder_name: String should match pattern '^[^/]+$'"],
        ),
        (
            "P1\tmale\tS1\tfalse\t\tL1\tWES\tF1\tHPO:1\t",
            ["2: hpo_terms: String should match pattern '^HP:\\d{7}(,HP:\\d{7})*$'"],
        ),
        (
            "P1\tmale\tS1\tfalse\t\tL1\tWES\tF1\t\t8M+T 8X+T",
            [
                "2: read_structure: String should match pattern '^((\\d+[TMBCS])*(\\d+|\\+)[TMBCS])( (\\d+[TMBCS])*(\\d+|\\+)[TMBCS])*$'"
            ],
        ),
        (
            "P1\tmale\tS1\tfalse\t\tL1\tRNA\tF1\t\t",
            ["2: library_type: 'RNA' is no known value"],
        ),
        ("P1\tmale\tS1\tfalse\t\tL1\t\tF1\t\t", ["2: library_type: a library needs its type"]),
        (
            "P1\tmale\tS1\tfalse\t\t\t\tF1\t\t",
            [
                "2: library_name: a row without a library describes a donor; it leaves sample_name, is_tumor, folder_name empty"
            ],
        ),
        (
            "P-1\tmale\tS1\tfalse\t\tL 1\tWES\tF1\t\t",
            ["2: library_name: String should match pattern '^[A-Za-z0-9_.-]+$'"],
        ),
    ],
)
def test_row_problems(tmp_path, row, problems):
    assert _problems(tmp_path, HEADER + row + "\n") == [f"samplesheet.tsv:{p}" for p in problems]


CORE = "family_id\tdonor_name\tfather\tmother\tsex\tis_index\tis_primary\tlibrary_name\tlibrary_type\tfolder_name\n"


def _germline(*rows):
    return CORE + "".join(row + "\n" for row in rows)


@pytest.mark.parametrize(
    "rows, problems",
    [
        (
            ["F\tP1\t\t\tmale\ttrue\t\tL1\tWES\tD1", "F\tP2\t\t\tmale\t\t\tL1\tWES\tD2"],
            ["3: library_name: 'L1' is also on line 2"],
        ),
        (
            ["F\tP1\t\t\tmale\ttrue\t\tL1\tWES\tD1", "F\tP2\t\t\tmale\t\t\tL2\tWES\tD1"],
            ["3: folder_name: 'D1' is also on line 2"],
        ),
        (
            ["F\tP1\t\t\tmale\t\t\tL1\tWES\tD1", "F\tP1\t\t\tfemale\t\t\tL2\tmRNA_seq\tD2"],
            ["3: sex: 'female' differs from 'male' on line 2 for donor 'P1'"],
        ),
        (
            ["F\tP1\t\t\tmale\t\t\tL1\tWES\tD1", "F\tP1\t\t\tmale\t\t\t\t\t"],
            [
                "3: library_name: donor 'P1' has a row without a library and other rows; a donor has library rows or one row without"
            ],
        ),
        (
            ["F\tP1\t0\t\tmale\t\t\tL1\tWES\tD1"],
            [
                "2: father: '0' is no donor of the sheets; a parent needs a row, with or without a library (0 is no value; leave the cell empty)"
            ],
        ),
        (["F\tP1\tP1\t\tmale\t\t\tL1\tWES\tD1"], ["2: father: a donor is not its own parent"]),
        (
            ["F\tP1\tP2\t\tmale\ttrue\t\tL1\tWES\tD1", "F\tP2\t\t\tfemale\t\t\tL2\tWES\tD2"],
            ["2: father: 'P2' is female (line 3)"],
        ),
        (
            ["F\tP1\t\tP2\tmale\t\t\tL1\tWES\tD1", "G\tP2\t\t\tfemale\t\t\tL2\tWES\tD2"],
            ["2: mother: 'P2' is in family 'G' (line 3), not in 'F'"],
        ),
        (
            [
                "F\tP1\tP3\t\tmale\ttrue\t\tL1\tWES\tD1",
                "F\tP2\t\tP3\tmale\t\t\tL2\tWES\tD2",
                "F\tP3\t\t\tunknown\t\t\tL3\tWES\tD3",
            ],
            ["4: donor_name: 'P3' is a father and a mother"],
        ),
        (
            ["F\tP1\tP2\t\tmale\ttrue\t\tL1\tWES\tD1", "F\tP2\tP1\t\tmale\t\t\tL2\tWES\tD2"],
            ["2: father: the parents of 'P1' form a cycle"],
        ),
        (
            ["F\tP1\t\t\tmale\t\t\tL1\tWES\tD1", "F\tP2\t\t\tmale\t\t\tL2\tWES\tD2"],
            [
                "2: is_index: family 'F' has 2 donors with a normal DNA library and no index; mark one is_index"
            ],
        ),
        (
            ["F\tP1\t\t\tmale\ttrue\t\tL1\tWES\tD1", "F\tP2\t\t\tmale\ttrue\t\tL2\tWES\tD2"],
            ["3: is_index: family 'F' has 2 index donors: P1, P2"],
        ),
        (
            ["F\tP1\t\t\tmale\ttrue\t\t\t\t", "F\tP2\t\t\tmale\t\t\tL2\tWES\tD2"],
            ["2: is_index: the index of a family has a normal DNA library"],
        ),
        (
            ["F\tP1\t\t\tmale\t\t\tL1\tWES\tD1", "F\tP1\t\t\tmale\t\t\tL2\tWGS\tD2"],
            [
                "3: is_primary: donor 'P1' has 2 normal dna libraries (L1, L2); mark exactly one is_primary, 0 are marked"
            ],
        ),
        (
            ["F\tP1\t\t\tmale\t\ttrue\tL1\tWES\tD1", "F\tP1\t\t\tmale\t\ttrue\tL2\tWGS\tD2"],
            [
                "3: is_primary: donor 'P1' has 2 normal dna libraries (L1, L2); mark exactly one is_primary, 2 are marked"
            ],
        ),
        (
            ["F\tP1\t\t\tmale\t\t\t\t\t"],
            ["samplesheet.tsv: library_name: the sample sheets have no library"],
        ),
    ],
)
def test_cross_row_problems(tmp_path, rows, problems):
    got = _problems(tmp_path, _germline(*rows))
    assert got == [p if p.startswith("samplesheet") else f"samplesheet.tsv:{p}" for p in problems]


def test_the_primary_library_is_the_marked_one_or_the_only_one(tmp_path):
    sheets = read_sheets(
        [
            _write(
                tmp_path,
                _germline(
                    "F\tP1\t\t\tmale\t\ttrue\tL1\tWES\tD1", "F\tP1\t\t\tmale\t\tfalse\tL2\tWGS\tD2"
                ),
            )
        ]
    )
    assert primary_library(sheets.libraries).values.library_name == "L1"
    assert primary_library(sheets.libraries[1:]).values.library_name == "L2"


def test_problems_are_sorted_and_capped(tmp_path):
    rows = "".join(f"P{i}\tmale\tS{i}\tfalse\t\tL{i}\tWES\tF{i}\tHPO\t\n" for i in range(60))
    with pytest.raises(SheetError) as excinfo:
        read_sheets([_write(tmp_path, HEADER + rows)])
    message = str(excinfo.value)
    assert message.startswith("The sample sheets have 60 problems:\n")
    assert message.count("\n  ") == 51 and message.endswith("... and 10 more")
    assert [p.line for p in excinfo.value.problems] == list(range(2, 62))


def test_columns_have_levels_and_isa_headers_and_the_header_lists_the_core_first():
    assert COLUMNS["donor_name"].level == Level.DONOR and COLUMNS["donor_name"].isa == "Source Name"
    assert (
        COLUMNS["is_tumor"].level == Level.SAMPLE and COLUMNS["library_kit"].level == Level.LIBRARY
    )
    assert header() == "donor_name\tlibrary_name\tlibrary_type\n"
    assert header(["is_tumor", "sex"]) == "donor_name\tsex\tis_tumor\tlibrary_name\tlibrary_type\n"


def test_synonyms_are_read_as_their_value_and_noted(tmp_path):
    text = (
        "donor_name\tsex\tis_affected\tis_tumor\tlibrary_name\tlibrary_type\tseq_platform\n"
        "P1\tF\tyes\tN\tL1\texome\tillumina\n"
        "P2\tXY\t\t1\tL2\tWGS\tPacBio\n"
    )
    sheets = read_sheets([_write(tmp_path, text)])
    rows = [row.values for row in sheets.rows]
    assert [(r.sex, r.is_affected, r.is_tumor, r.library_type, r.seq_platform) for r in rows] == [
        ("female", True, False, "WES", "Illumina"),
        ("male", None, True, "WGS", "PacBio"),
    ]
    assert [str(n).split(":", 1)[1] for n in sheets.notes] == [
        "2: sex: 'F' is read as 'female'",
        "2: is_affected: 'yes' is read as 'true'",
        "2: is_tumor: 'N' is read as 'false'",
        "2: library_type: 'exome' is read as 'WES'",
        "2: seq_platform: 'illumina' is read as 'Illumina'",
        "3: sex: 'XY' is read as 'male'",
        "3: is_tumor: '1' is read as 'true'",
    ]


def test_ped_codes_of_is_affected_are_rejected_and_unknown_values_too(tmp_path):
    text = "donor_name\tsex\tis_affected\tlibrary_name\tlibrary_type\nP1\tdiverse\t2\tL1\tWES\n"
    assert _problems(tmp_path, text) == [
        "samplesheet.tsv:2: is_affected: '2' is no known value; the PED codes 1 and 2 are ambiguous, write true or false",
        "samplesheet.tsv:2: sex: 'diverse' is no known value",
    ]


def test_fix_sheet_rewrites_the_synonyms_and_nothing_else(tmp_path):
    from snappy_pipeline.sample_sheet import fix_sheet

    text = (
        "# kept\r\n"
        "donor_name\tsex\tis_tumor\tlibrary_name\tlibrary_type\tnote\r\n"
        "\r\n"
        "P1\tF\tY\tL1\twes\tkeep F\r\n"
        "P2\tmale\tfalse\tL2\tWGS\t\r\n"
    )
    path = _write(tmp_path, text)
    assert fix_sheet(path) == 3
    assert (tmp_path / "samplesheet.tsv").read_bytes().decode() == (
        "# kept\r\n"
        "donor_name\tsex\tis_tumor\tlibrary_name\tlibrary_type\tnote\r\n"
        "\r\n"
        "P1\tfemale\ttrue\tL1\tWES\tkeep F\r\n"
        "P2\tmale\tfalse\tL2\tWGS\t\r\n"
    )
    assert fix_sheet(path) == 0 and read_sheets([path]).notes == []


def test_the_platform_follows_from_the_instrument(tmp_path):
    text = (
        "donor_name\tlibrary_name\tlibrary_type\tinstrument_model\tseq_platform\n"
        "P1\tL1\tWES\tNextSeq 1000\t\n"
        "P2\tL2\tWGS\tRevio\t\n"
        "P3\tL3\tWGS\tPromethION 2\t\n"
        "P4\tL4\tWGS\tAVITI\t\n"
        "P5\tL5\tWGS\tUnheard-of 3000\t\n"
        "P6\tL6\tWGS\tNovaSeq X\tother\n"
    )
    rows = [row.values for row in read_sheets([_write(tmp_path, text)]).rows]
    assert [row.seq_platform for row in rows] == [
        "Illumina",
        "PacBio",
        "ONT",
        "other",
        None,
        "other",
    ]
    assert rows[0].instrument_model == "NextSeq 1000"


def test_a_sequencer_in_seq_platform_names_the_instrument_column(tmp_path):
    text = "donor_name\tlibrary_name\tlibrary_type\tseq_platform\nP1\tL1\tWES\tNextSeq 1000\n"
    assert _problems(tmp_path, text) == [
        "samplesheet.tsv:2: seq_platform: 'NextSeq 1000' is no known value; the platform is "
        "Illumina, PacBio, ONT or other, a sequencer such as 'NextSeq 1000' goes into "
        "instrument_model, which gives the platform"
    ]
