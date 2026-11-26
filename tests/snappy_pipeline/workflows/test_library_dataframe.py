# -*- coding: utf-8 -*-
"""Tests for the library table built from the plain sample sheets."""

from pathlib import Path

import pandas as pd
import pytest

from snappy_pipeline.base import ConfigError
from snappy_pipeline.library_table import (
    build_donor_dataframe,
    build_library_dataframe,
    entity_libraries,
    resolve_related,
    resolve_relationships,
    select_libraries,
)
from snappy_pipeline.sample_sheet import read_sheets

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _table(tmp_path, text, include_background=False, name="s.tsv"):
    (tmp_path / name).write_text(text)
    return build_library_dataframe(
        read_sheets([str(tmp_path / name)]), include_background=include_background
    )


def _sheets(tmp_path, text, name="s.tsv"):
    (tmp_path / name).write_text(text)
    return read_sheets([str(tmp_path / name)])


def test_cancer_sheet_rows(tmp_path):
    df = _table(tmp_path, (FIXTURES / "samplesheet_rna.tsv").read_text())
    assert df["library_name"].tolist() == [
        "case001subregion-N1-DNA1-WES1",
        "case001subregion-T1-DNA1-WES1",
        "case001subregion-T1-RNA1-mRNA_seq1",
    ]
    assert df["role"].tolist() == ["normal", "tumor", "tumor"]
    assert df["extraction_type"].tolist() == ["dna", "dna", "rna"]
    assert df["is_primary"].tolist() == [True, True, True]  # alone in their role and type
    assert df["is_index"].tolist() == [True] * 3  # the only donor of its family
    assert df["sex"].tolist() == ["unknown"] * 3
    assert df["father_name"].tolist() == ["0"] * 3 and df["disease_state"].tolist() == [0] * 3
    assert pd.isna(df["library_kit"].tolist()[2]) and df["is_background"].tolist() == [False] * 3
    assert "kind" not in df.columns and "related_group" not in df.columns


def test_germline_trio_rows(germline_sheet_tsv, tmp_path):
    df = _table(tmp_path, germline_sheet_tsv)
    p001 = df.set_index("donor_name").loc["P001"]
    assert p001["family_id"] == "P001" and p001["is_index"] and p001["disease_state"] == 2
    assert p001["father_name"] == "P002-N1-DNA1-WGS1" and p001["mother_name"] == "P003-N1-DNA1-WGS1"
    assert df.set_index("donor_name").loc["P002"]["disease_state"] == 1
    assert df["is_index"].tolist() == [True, False, False, True, False, False]
    assert df["sex"].tolist() == ["female", "male", "female", "male", "male", "female"]


def test_a_donors_other_libraries_and_the_primary_one(tmp_path):
    text = (
        "donor_name\tsex\tlibrary_name\tlibrary_type\tis_primary\n"
        "P1\tmale\tP1-WES\tWES\tfalse\n"
        "P1\tmale\tP1-WGS\tWGS\ttrue\n"
        "P1\tmale\tP1-RNA\tmRNA_seq\t\n"
        "P2\tfemale\tP2-WES\tWES\t\n"
    )
    df = _table(tmp_path, text)
    assert df["is_primary"].tolist() == [False, True, True, True]
    donors = build_donor_dataframe(read_sheets([str(tmp_path / "s.tsv")]))
    assert donors["library_name"].tolist() == ["P1-WGS", "P2-WES"]
    assert donors["pedigree"].tolist() == ["", ""]


def test_relapse_sheet_keeps_the_extra_column_and_the_primary_tumor(tmp_path):
    df = _table(tmp_path, (FIXTURES / "samplesheet_somatic_trio.tsv").read_text())
    assert pd.isna(df["timepoint"][0]) and df["timepoint"].tolist()[1:] == ["primary", "relapse"]
    assert df["is_primary"].tolist() == [True, True, False]
    df = resolve_relationships(
        df, {"matched_normal_lib": {"via": "donor_name", "target": "role == 'normal'"}}
    )
    assert df["matched_normal_lib"].tolist() == ["case001subregion-N1-DNA1-WES1"] * 3


def test_background_libraries_are_left_out_unless_asked_for(tmp_path):
    text = (
        "donor_name\tlibrary_name\tlibrary_type\tis_background\n"
        "P1\tL1\tWES\t\n"
        "C1\tC1-L\tWES\ttrue\n"
    )
    assert _table(tmp_path, text)["library_name"].tolist() == ["L1"]
    df = _table(tmp_path, text, include_background=True)
    assert df["is_background"].tolist() == [False, True]
    donors = build_donor_dataframe(_sheets(tmp_path, text), include_background=False)
    assert donors["donor_name"].tolist() == ["P1"]


def test_the_donor_table_lists_unsequenced_donors_with_their_parents(tmp_path):
    text = (
        "family_id\tdonor_name\tfather\tmother\tsex\tis_index\tlibrary_name\tlibrary_type\n"
        "F\tidx\tfa\tmo\tfemale\ttrue\tidx-L\tWES\n"
        "F\tfa\t\t\tmale\t\tfa-L\tWES\n"
        "F\tmo\t\t\tfemale\t\t\t\n"
    )
    donors = build_donor_dataframe(_sheets(tmp_path, text))
    assert donors.to_dict("records") == [
        {
            "donor_name": "idx",
            "father_name": "fa",
            "mother_name": "mo",
            "sex": "female",
            "library_name": "idx-L",
            "pedigree": "F",
            "is_background": False,
        },
        {
            "donor_name": "fa",
            "father_name": "0",
            "mother_name": "0",
            "sex": "male",
            "library_name": "fa-L",
            "pedigree": "F",
            "is_background": False,
        },
        {
            "donor_name": "mo",
            "father_name": "0",
            "mother_name": "0",
            "sex": "female",
            "library_name": "",
            "pedigree": "F",
            "is_background": False,
        },
    ]
    df = _table(tmp_path, text)
    assert df.set_index("donor_name").loc["idx"]["mother_name"] == "0"  # no library


def test_default_selection_and_requirements(tmp_path):
    df = _table(tmp_path, (FIXTURES / "samplesheet_rna.tsv").read_text())
    assert select_libraries(df, None)["extraction_type"].tolist() == ["dna", "dna"]
    assert select_libraries(df, None, requirement="role == 'tumor'")["library_name"].tolist() == [
        "case001subregion-T1-DNA1-WES1"
    ]
    assert select_libraries(df, "extraction_type == 'rna'", default="")[
        "library_name"
    ].tolist() == ["case001subregion-T1-RNA1-mRNA_seq1"]
    assert len(select_libraries(df, None, default="")) == 3
    with pytest.raises(ConfigError, match="library_selection expression 'nope ==' is invalid"):
        select_libraries(df, "nope ==")


def test_entities_by_library_donor_family_and_related(germline_trio_plus_sheet_tsv, tmp_path):
    df = _table(tmp_path, germline_trio_plus_sheet_tsv)
    assert len(entity_libraries(df, "library")) == 7
    assert entity_libraries(df, "family") == {
        "family1": "P001-N1-DNA1-WGS1",
        "family2": "P004-N1-DNA1-WGS1",
    }
    assert list(entity_libraries(df, "donor")) == df["donor_name"].tolist()
    assert resolve_related(df, "related") == "family"
    cancer = _table(tmp_path, (FIXTURES / "samplesheet_somatic_trio.tsv").read_text(), name="c.tsv")
    assert resolve_related(cancer, "related") == "donor"
    assert entity_libraries(cancer, "related") == {
        "case001subregion": "case001subregion-N1-DNA1-WES1"
    }
    mixed = _table(
        tmp_path,
        "family_id\tdonor_name\tlibrary_name\tlibrary_type\nF\tP1\tL1\tWES\n\tP2\tL2\tWES\n",
        name="m.tsv",
    )
    with pytest.raises(
        ConfigError,
        match="group_by: related is the family .* set group_by: family or group_by: donor",
    ):
        resolve_related(mixed, "related")
    with pytest.raises(ConfigError, match="group_by: family needs a family_id .* L2"):
        entity_libraries(mixed, "family")
