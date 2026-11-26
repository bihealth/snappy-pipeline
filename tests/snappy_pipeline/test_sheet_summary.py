# -*- coding: utf-8 -*-
"""Tests for the check, summary and warnings of the libraries in the sample sheets."""

from types import SimpleNamespace

import pandas as pd
import pytest
from snappy_pipeline.sheet_summary import (
    check_libraries,
    library_summary,
    library_warnings,
    missing_library_files,
)


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


def test_summary_splits_a_data_set_by_the_batches_of_the_sheet():
    df = _libraries(
        ("a", "first", "germline", "dna", "WES", "MedExome"),
        ("b", "first", "germline", "dna", "WES", "MedExome"),
        ("c", "first", "germline", "dna", "WES", "MedExome"),
    )
    df["batch_no"] = ["run-2024-A", "run-2024-B", "run-2024-A"]
    lines = library_summary(df).splitlines()
    assert lines[0].split()[:3] == ["data_set", "batch_no", "kind"]
    assert lines[1].split() == ["first", "run-2024-A", "germline", "dna", "WES", "MedExome", "2"]
    assert lines[2].split() == ["first", "run-2024-B", "germline", "dna", "WES", "MedExome", "1"]


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


def _project(tmp_path, table_libraries=("a", "b")):
    """A project with the external_data tasks raw (reads) and tables (BAMs), whose files tables
    exist and list ``table_libraries``, and missing, whose table does not exist."""
    steps = {"raw": "external_data", "tables": "external_data", "missing": "external_data"}
    (tmp_path / "r.tsv").touch()
    (tmp_path / "t.tsv").touch()
    return SimpleNamespace(
        task=lambda name: SimpleNamespace(step=steps.get(name, "trimming")),
        dependencies={},
        task_configs={
            "raw": SimpleNamespace(
                library_files=[str(tmp_path / "r.tsv")], patterns_find_reads=lambda: True
            ),
            "tables": SimpleNamespace(
                library_files=[str(tmp_path / "t.tsv")], patterns_find_reads=lambda: False
            ),
            "missing": SimpleNamespace(
                library_files=["/none/m.tsv"], patterns_find_reads=lambda: True
            ),
        },
        files_table=lambda config: SimpleNamespace(
            files=lambda library: {"bam": "a.bam"} if library == "a" else {},
            libraries=list(table_libraries),
        ),
    )


def _task(project, name, field, source, libraries, read_groups=None):
    project.dependencies[name] = {field: source}
    return SimpleNamespace(
        depends_on=SimpleNamespace(reads=source if field == "reads" else ""),
        project=project,
        libraries=pd.DataFrame({"library_name": libraries}),
        read_groups=lambda library, source=None: read_groups(library),
        reads_from_task=lambda source: project.task(source).step != "external_data",
        library_table=lambda include_background: pd.DataFrame({"library_name": ["a", "b"]}),
    )


def test_missing_reads_checks_each_library_of_a_source_once_and_skips_trimmed_reads(tmp_path):
    def read_groups(library):
        if library == "b":
            raise ValueError("The files tables r.tsv list no files of library 'b'")
        return ["group"]

    def fail(library):
        raise AssertionError("the trimming task finds the reads, not its reader")

    project = _project(tmp_path)
    tasks = {
        "mapping": _task(project, "mapping", "reads", "raw", ["a", "b"], read_groups),
        "qc": _task(project, "qc", "reads", "raw", ["b"], read_groups),
        "trimmed": _task(project, "trimmed", "reads", "trimming", ["b"], fail),
        "genome": SimpleNamespace(depends_on=None, project=project),
    }
    assert missing_library_files(tasks) == [
        ("mapping", "b", "The files tables r.tsv list no files of library 'b'"),
    ]


def test_missing_library_files_checks_files_tables_and_their_unknown_libraries(tmp_path):
    project = _project(tmp_path, table_libraries=("a", "zzz"))
    tasks = {"calling": _task(project, "calling", "alignments", "tables", ["a", "b"])}
    table = tmp_path / "t.tsv"
    assert missing_library_files(tasks) == [
        (
            "calling",
            "b",
            f"The files tables {table} list no files of library 'b'; add them to the table, or "
            "run `snappy scan-files --task tables --pattern key=regex ...`",
        ),
        ("tables", "zzz", "is in a files table but in no sample sheet"),
    ]


def test_a_missing_files_table_is_reported_once_with_the_command_that_writes_it(tmp_path):
    project = _project(tmp_path)
    tasks = {
        "mapping": _task(project, "mapping", "reads", "missing", ["a", "b"]),
        "qc": _task(project, "qc", "reads", "missing", ["a"]),
    }
    assert missing_library_files(tasks) == [
        ("missing", "/none/m.tsv", "does not exist; run `snappy scan-files`"),
    ]
