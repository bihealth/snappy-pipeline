# -*- coding: utf-8 -*-
"""Check of the libraries in the sample sheets, with a summary and warnings.

The main Snakefile runs :func:`print_library_summary` once per run, in the main Snakemake
process, so cluster jobs do not repeat it. It ends with :func:`missing_library_files`, which
``snappy check-files`` runs as well: every library of a task that reads FASTQ files must have
reads, or the DAG build would fail later with the first of them.
"""

import os
from collections.abc import Mapping
from typing import Any

import pandas as pd
from jsonschema import Draft202012Validator
from ruamel.yaml import YAML

from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.external_data import no_rows_message
from snappy_pipeline.library_table import TARGETED_LIBRARY_TYPES

#: JSON schema of the sample sheet values of one library
LIBRARY_SCHEMA = os.path.join(os.path.dirname(__file__), "schemas", "library.schema.yaml")

#: Columns of the summary table, where the sample sheets have them (``batch_no`` is the sheet's
#: ``batchNo`` custom field)
SUMMARY_COLUMNS = ["data_set", "batch_no", "kind", "extraction_type", "library_type", "library_kit"]

#: Number of libraries that a warning names before it shortens the list
MAX_NAMES = 5


def _has_value(value: Any) -> bool:
    return value is not None and value == value  # NaN is the only value unequal to itself


def check_libraries(df: pd.DataFrame) -> None:
    """Raise :class:`ConfigError` naming every library whose values break the schema.

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
        raise ConfigError("Sample sheet values out of range:\n" + "\n".join(errors))


def library_summary(df: pd.DataFrame) -> str:
    """Return the number of libraries per data set, batch, kind, extraction type, library type and kit."""
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


def missing_library_files(instances: Mapping[str, BaseStep]) -> list[tuple[str, str, str]]:
    """Return ``(task, library, problem)`` for every library of a task whose files an
    ``external_data`` source does not provide: the files tables of the task (``library_files``)
    have no row for the library. Rows of such a table for libraries that no sample sheet has are
    reported for the table's task.

    A task that reads the FASTQs of another task (trimmed reads) is covered by that task. Each
    library is checked once per source, with the first task that reads it.
    """
    problems = []
    checked: set[tuple[str, str]] = set()
    table_tasks: dict[str, BaseStep] = {}
    for name, task in instances.items():
        for field, source in task.project.dependencies.get(name, {}).items():
            if not isinstance(source, str):
                continue
            config = task.project.task_configs.get(source)
            reads = field == "reads"
            tables = getattr(config, "library_files", None)
            if reads and task.reads_from_task(source):
                continue
            if not reads and not tables:
                continue
            if tables:
                if missing := [table for table in tables if not os.path.exists(table)]:
                    fix = "run `snappy scan-files`" if config.scan is not None else "write it"
                    problems += [
                        (source, table, f"does not exist; {fix}")
                        for table in missing
                        if (source, table) not in checked
                    ]
                    checked.update((source, table) for table in missing)
                    continue
                table_tasks.setdefault(source, task)
            for library in task.libraries["library_name"]:
                if (source, library) in checked:
                    continue
                checked.add((source, library))
                try:
                    if reads:
                        task.read_groups(library, source)
                    elif not task.project.files_table(config).files(library):
                        problems.append((name, library, no_rows_message(config, library)))
                except (ValueError, ConfigError) as e:
                    problems.append((name, library, str(e)))
    for source, task in table_tasks.items():
        config = task.project.task_configs[source]
        known = set(task.library_table(include_background=True)["library_name"])
        table = task.project.files_table(config)
        for library in table.libraries:
            if library not in known:
                problems.append((source, library, "is in a files table but in no sample sheet"))
    return problems


def print_library_summary(instances: Mapping[str, BaseStep]) -> None:
    """Check the libraries of the sample sheets, print their summary and the warnings of the
    sample sheets and of the tasks, then raise :class:`ConfigError` if a library has no reads."""
    if not instances:
        return
    df = next(iter(instances.values())).library_table(include_background=True)
    check_libraries(df)
    print(f"[snappy] Libraries in the sample sheets:\n{library_summary(df)}", flush=True)
    warnings = library_warnings(df)
    warnings += [warning for task in instances.values() for warning in task.library_warnings()]
    for warning in warnings:
        print(f"[snappy] Warning: {warning}", flush=True)
    if problems := missing_library_files(instances):
        lines = [f"{task}: {library}: {problem}" for task, library, problem in problems]
        raise ConfigError("Libraries without files:\n  " + "\n  ".join(lines))
