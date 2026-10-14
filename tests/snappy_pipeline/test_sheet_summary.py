# -*- coding: utf-8 -*-
"""Tests for the check, summary and warnings of the libraries in the sample sheets."""

import pandas as pd
import pytest
from snappy_pipeline.sheet_summary import check_libraries, library_summary, library_warnings


def _libraries(*rows):
    columns = ["library_name", "data_set", "kind", "extraction_type", "library_type", "library_kit"]
    return pd.DataFrame([dict(zip(columns, row)) for row in rows], columns=columns)


def test_summary_counts_libraries_per_data_set_type_and_kit():
    df = _libraries(
        ("a", "first", "germline", "dna", "WES", "MedExome"),
        ("b", "first", "germline", "dna", "WES", "MedExome"),
        ("c", "second", "cancer", "rna", "mRNA_seq", None),
    )
    lines = library_summary(df).splitlines()
    assert lines[0].split() == [
        "data_set",
        "kind",
        "extraction_type",
        "library_type",
        "library_kit",
        "libraries",
    ]
    assert lines[1].split() == ["first", "germline", "dna", "WES", "MedExome", "2"]
    assert lines[2].split() == ["second", "cancer", "rna", "mRNA_seq", "1"]


def test_warnings_name_mixed_rna_types_several_kits_and_missing_kits():
    df = _libraries(
        ("a", "first", "germline", "dna", "WES", "MedExome"),
        ("b", "first", "germline", "dna", "WES", "Twist_Exome_2"),
        ("c", "first", "germline", "dna", "Panel_seq", None),
        ("d", "first", "germline", "dna", "WGS", "TruSeq_PCR_free"),
        ("e", "first", "germline", "rna", "mRNA_seq", None),
        ("f", "first", "germline", "rna", "total_RNA_seq", None),
    )
    assert library_warnings(df) == [
        "Data set 'first' has mRNA_seq and total_RNA_seq libraries",
        "Data set 'first' has targeted libraries of several kits: MedExome (1), Twist_Exome_2 (1)",
        "Data set 'first' has 1 targeted libraries without libraryKit: c",
    ]


def test_no_warnings_for_one_kit():
    df = _libraries(
        ("a", "first", "germline", "dna", "WES", "MedExome"),
        ("b", "first", "germline", "dna", "WES", "MedExome"),
    )
    assert library_warnings(df) == []


def test_schema_accepts_long_read_chemistry_as_library_kit():
    check_libraries(_libraries(("a", "first", "germline", "dna", "WGS", "PacBio HiFi")))


@pytest.mark.parametrize(
    "row",
    [
        ("a", "first", "germline", "unknown", "WES", "MedExome"),
        ("a", "first", "germline", "dna", "WES", "SureSelect V6"),
        ("a", "first", "germline", "dna", "Panel_seq", "kit/1"),
    ],
)
def test_schema_rejects_unknown_extraction_type_and_bad_kit_names(row):
    with pytest.raises(ValueError, match="'a' of data set 'first'"):
        check_libraries(_libraries(row))
