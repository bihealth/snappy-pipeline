# -*- coding: utf-8 -*-
"""Check of the libraries in the sample sheets, with a summary and warnings.

The main Snakefile runs :func:`print_library_summary` once per run, in the main Snakemake
process, so cluster jobs do not repeat it.
"""

import os
from collections.abc import Mapping
from typing import Any

import pandas as pd
from jsonschema import Draft202012Validator
from ruamel.yaml import YAML

from snappy_pipeline.workflows.abstract import TARGETED_LIBRARY_TYPES, BaseStep

#: JSON schema of the sample sheet values of one library
LIBRARY_SCHEMA = os.path.join(os.path.dirname(__file__), "schemas", "library.schema.yaml")

#: Columns of the summary table, where the sample sheets have them
SUMMARY_COLUMNS = ["data_set", "kind", "extraction_type", "library_type", "library_kit"]

#: Number of libraries that a warning names before it shortens the list
MAX_NAMES = 5


def _has_value(value: Any) -> bool:
    return value is not None and value == value  # NaN is the only value unequal to itself


def check_libraries(df: pd.DataFrame) -> None:
    """Raise ``ValueError`` naming every library whose values break the schema.

    Not ``snakemake.utils.validate``: it needs jsonschema 4.18 or newer.
    """
    with open(LIBRARY_SCHEMA) as f:
        validator = Draft202012Validator(YAML(typ="safe").load(f))
    errors = []
    for record in df.to_dict("records"):
        record = {key: value for key, value in record.items() if _has_value(value)}
        name = f"{record.get('library_name')!r} of data set {record.get('data_set')!r}"
        errors += [f"{name}: {error.message}" for error in validator.iter_errors(record)]
    if errors:
        raise ValueError("Sample sheet values out of range:\n" + "\n".join(errors))


def library_summary(df: pd.DataFrame) -> str:
    """Return the number of libraries per data set, kind, extraction type, library type and kit."""
    columns = [column for column in SUMMARY_COLUMNS if column in df.columns]
    counts = df[columns].fillna("").groupby(columns).size().rename("libraries")
    return counts.reset_index().to_string(index=False)


def _names(names: list[str]) -> str:
    shown = ", ".join(names[:MAX_NAMES])
    return shown + (", ..." if len(names) > MAX_NAMES else "")


def library_warnings(df: pd.DataFrame) -> list[str]:
    """Return warnings about the libraries of each data set: mRNA and total RNA libraries
    together, targeted libraries of several kits, and targeted libraries without a kit."""
    warnings = []
    if "library_type" not in df.columns:
        return warnings
    for data_set, libraries in df.groupby("data_set"):
        if {"mRNA_seq", "total_RNA_seq"} <= set(libraries["library_type"]):
            warnings.append(f"Data set {data_set!r} has mRNA_seq and total_RNA_seq libraries")
        targeted = libraries[libraries["library_type"].isin(TARGETED_LIBRARY_TYPES)]
        kits = targeted.get("library_kit", pd.Series("", index=targeted.index)).fillna("")
        counts = kits[kits != ""].value_counts()
        if len(counts) > 1:
            listed = ", ".join(f"{kit} ({n})" for kit, n in sorted(counts.items()))
            warnings.append(
                f"Data set {data_set!r} has targeted libraries of several kits: {listed}"
            )
        if without := targeted["library_name"][kits == ""].tolist():
            warnings.append(
                f"Data set {data_set!r} has {len(without)} targeted libraries without libraryKit: "
                + _names(without)
            )
    return warnings


def print_library_summary(instances: Mapping[str, BaseStep]) -> None:
    """Check the libraries of the sample sheets, then print their summary and the warnings of
    the sample sheets and of the tasks."""
    if not instances:
        return
    df = next(iter(instances.values())).build_library_dataframe(include_background=True)
    check_libraries(df)
    print(f"[snappy] Libraries in the sample sheets:\n{library_summary(df)}", flush=True)
    warnings = library_warnings(df)
    warnings += [warning for task in instances.values() for warning in task.library_warnings()]
    for warning in warnings:
        print(f"[snappy] Warning: {warning}", flush=True)
