# -*- coding: utf-8 -*-
"""Check of the libraries in the sample sheets, with a summary and warnings.

The main Snakefile runs :func:`print_library_summary` once per run, in the main Snakemake
process, so cluster jobs do not repeat it. It ends with :func:`missing_library_files`, which
``snappy check-files`` runs as well: every library of a task that reads FASTQ files must have
reads, or the DAG build would fail later with the first of them.
"""

import os
from collections.abc import Mapping

import pandas as pd

from snappy_pipeline.base import ConfigError
from snappy_pipeline.library_table import TARGETED_LIBRARY_TYPES
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.external_data import no_rows_message

#: Columns of the summary table, where the sample sheets have values
SUMMARY_COLUMNS = ["batch_no", "role", "extraction_type", "library_type", "library_kit"]

#: Number of libraries that a warning names before it shortens the list
MAX_NAMES = 5


def library_summary(df: pd.DataFrame) -> str:
    """Return the number of libraries per batch, role, extraction type, library type and kit."""
    columns = [c for c in SUMMARY_COLUMNS if c in df.columns and df[c].notna().any()]
    counts = df[columns].fillna("").groupby(columns).size().rename("libraries")
    return counts.reset_index().to_string(index=False)


def _names(names: list[str]) -> str:
    shown = ", ".join(names[:MAX_NAMES])
    return shown + (", ..." if len(names) > MAX_NAMES else "")


def library_warnings(df: pd.DataFrame) -> list[str]:
    """Return warnings about the libraries: mRNA and total RNA libraries together, targeted
    libraries of several kits, and targeted libraries without a kit."""
    warnings = []
    if df.empty:
        return warnings
    if {"mRNA_seq", "total_RNA_seq"} <= set(df["library_type"]):
        warnings.append("The sample sheets have mRNA_seq and total_RNA_seq libraries")
    targeted = df[df["library_type"].isin(TARGETED_LIBRARY_TYPES)]
    kits = targeted["library_kit"].fillna("")
    counts = kits[kits != ""].value_counts()
    if len(counts) > 1:
        listed = ", ".join(f"{kit} ({n})" for kit, n in sorted(counts.items()))
        warnings.append(f"The targeted libraries have several kits: {listed}")
    if without := targeted["library_name"][kits == ""].tolist():
        warnings.append(
            f"{len(without)} targeted libraries have no library_kit: " + _names(without)
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
                    fix = (
                        "run `snappy scan-files`"
                        if config.provides_reads()
                        else f"run `snappy scan-files --task {source} --pattern key=regex ...`"
                    )
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
                        problems.append((name, library, no_rows_message(config, library, source)))
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
    print(f"[snappy] Libraries in the sample sheets:\n{library_summary(df)}", flush=True)
    warnings = library_warnings(df)
    warnings += [warning for task in instances.values() for warning in task.library_warnings()]
    for warning in warnings:
        print(f"[snappy] Warning: {warning}", flush=True)
    if problems := missing_library_files(instances):
        lines = [f"{task}: {library}: {problem}" for task, library, problem in problems]
        raise ConfigError("Libraries without files:\n  " + "\n  ".join(lines))
