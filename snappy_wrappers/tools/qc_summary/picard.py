# -*- coding: utf-8 -*-
"""The Picard metrics files of an ``ngs_data_qc`` task with ``tool: picard``.

The files come by program (the kinds of ``files``). Every ``METRICS CLASS`` section becomes a
table with one row per library (and per ``CATEGORY`` or ``PAIR_ORIENTATION`` where Picard
reports several); every known ``HISTOGRAM`` section becomes a chart (:data:`HISTOGRAMS`). The
summary table joins the programs' metrics into one row per library, the columns prefixed with
the metrics class.
"""

from __future__ import annotations

from functools import reduce

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    Chart,
    Column,
    Summary,
    Table,
    metric_strips,
    numeric,
    thin,
    ticks,
)

#: Charts of the histogram sections, by their first column: (name, title, y column or ``None``
#: for the sum of the other columns, y title, log y)
HISTOGRAMS = {
    "coverage_or_base_quality": (
        "target_coverage",
        "Target coverage distribution",
        "high_quality_coverage_count",
        "target bases",
        True,
    ),
    "coverage": ("coverage", "Coverage distribution", "high_quality_coverage_count", "bases", True),
    "insert_size": ("insert_size", "Insert size distribution", None, "read pairs", False),
    "CYCLE": ("base_quality", "Mean base quality by cycle", "MEAN_QUALITY", "mean quality", False),
    "QUALITY": ("quality_distribution", "Base quality distribution", "COUNT_OF_Q", "bases", False),
    "normalized_position": (
        "gene_body_coverage",
        "Gene body coverage",
        "All_Reads.normalized_coverage",
        "normalized coverage",
        False,
    ),
}

#: The row of a library that the summary table keeps when a program reports several
PREFERRED_ROWS = {"CATEGORY": "PAIR", "PAIR_ORIENTATION": "FR", "ACCUMULATION_LEVEL": "ALL READS"}

#: Text columns of the metrics tables
LABEL_COLUMNS = ("BAIT_SET", "SAMPLE", "LIBRARY", "READ_GROUP", *PREFERRED_ROWS)


def parse_metrics(path: str) -> tuple[str | None, pl.DataFrame | None, pl.DataFrame | None]:
    """Return the metrics class, the metrics table and the histogram of a Picard metrics file."""
    metrics_class, metrics, histogram = None, None, None
    section: list[str] = []
    kind = None
    with open(path, "rt") as f:
        lines = f.read().splitlines()
    for line in [*lines, ""]:
        if line.startswith("## METRICS CLASS"):
            kind, section = "metrics", []
            metrics_class = line.split("\t", 1)[1].rsplit(".", 1)[-1]
            continue
        if line.startswith("## HISTOGRAM"):
            kind, section = "histogram", []
            continue
        if line.startswith("#"):
            continue
        if line.strip():
            if kind:
                section.append(line)
            continue
        if kind and len(section) > 1:
            if kind == "metrics":
                metrics = _table(section)
            else:
                histogram = _table(section)
        kind, section = None, []
    return metrics_class, metrics, histogram


def _table(lines: list[str]) -> pl.DataFrame:
    """A frame of the tab-separated ``lines``, the first the header; numbers where possible."""
    header = lines[0].split("\t")
    rows = [line.split("\t")[: len(header)] for line in lines[1:]]
    frame = pl.DataFrame(rows, schema=header, orient="row")
    return numeric(frame, [c for c in header if c not in LABEL_COLUMNS])


def summarize(files: dict[str, list[str]], labels: list[str]) -> Summary:
    tables: dict[str, list[pl.DataFrame]] = {}
    histograms: dict[str, list[pl.DataFrame]] = {}
    for program, paths in files.items():
        for label, path in zip(labels, paths):
            metrics_class, metrics, histogram = parse_metrics(path)
            library = pl.lit(label).alias("library")
            if metrics is not None:
                tables.setdefault(metrics_class or program, []).append(
                    metrics.with_columns(library)
                )
            if histogram is not None and histogram.columns[0] in HISTOGRAMS:
                histograms.setdefault(histogram.columns[0], []).append(
                    histogram.with_columns(library)
                )

    views: list[Table | Chart] = []
    strips: list[Table | Chart] = []
    summary_parts = []
    for name, frames in tables.items():
        frame = pl.concat(frames, how="diagonal")
        frame = frame.select("library", pl.exclude("library"))
        views.append(
            Table(
                name,
                name,
                frame,
                desc=f"Picard `{name}` of each library.",
                columns={c: _column(c) for c in frame.columns if c != "library"},
            )
        )
        part = _one_row_per_library(frame)
        summary_parts.append(
            part.rename({c: f"{name}.{c}" for c in part.columns if c != "library"})
        )
    for first_column, frames in histograms.items():
        name, title, y, y_title, log_y = HISTOGRAMS[first_column]
        frame = pl.concat(frames, how="diagonal")
        x = first_column
        if y is None:
            counts = [c for c in frame.columns if c not in ("library", x)]
            frame = frame.with_columns(pl.sum_horizontal(counts).alias(y_title))
            y = y_title
        elif y not in frame.columns:
            continue
        frame = frame.select("library", x, y).drop_nulls()
        views.append(
            Chart(
                name,
                title,
                thin(frame, x, y),
                x=x,
                y=y,
                x_title=x.replace("_", " ").lower(),
                y_title=y_title,
                log_y=log_y,
            )
        )

    if summary_parts:
        summary = reduce(
            lambda a, b: a.join(b, on="library", how="full", coalesce=True), summary_parts
        )
    else:
        summary = pl.DataFrame({"library": labels})
    columns = {c: _column(c.rsplit(".", 1)[-1]) for c in summary.columns if c != "library"}
    metrics = {
        c: c.rsplit(".", 1)[-1].replace("_", " ").lower()
        for c in summary.columns
        if c != "library" and _column(c.rsplit(".", 1)[-1]).plot is not None
    }
    if metrics:
        strips.append(
            metric_strips(
                summary,
                metrics,
                "libraries",
                "Libraries by metric",
                desc="One tick per library on the range of each fraction, mean and median: a "
                "tick apart from the others is a library to look at; the tooltip names it.",
            )
        )
    return Summary(
        table=Table(
            "summary",
            "Picard metrics",
            summary,
            desc="The Picard metrics of each library, one row per library, columns prefixed with "
            "the metrics class.",
            columns=columns,
        ),
        views=[*strips, *views],
    )


def _one_row_per_library(frame: pl.DataFrame) -> pl.DataFrame:
    """Keep one row per library: the preferred category, orientation or level, else the first."""
    for column, value in PREFERRED_ROWS.items():
        if column in frame.columns and (frame[column] == value).any():
            frame = frame.filter(pl.col(column) == value)
            break
    return frame.unique(subset=["library"], keep="first", maintain_order=True)


def _column(name: str) -> Column:
    if name.startswith("PCT_") or name.endswith(("_FRACTION", "_RATE", "_PCT")):
        return ticks()
    if name.startswith(("MEAN_", "MEDIAN_")):
        return ticks(precision=1)
    if name in LABEL_COLUMNS:
        return Column()
    return Column(detail=True)
