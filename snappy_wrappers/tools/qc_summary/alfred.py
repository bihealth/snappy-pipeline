# -*- coding: utf-8 -*-
"""The ``alfred qc`` JSON files of an ``ngs_data_qc`` task with ``tool: alfred``.

The summary table has one row per library and read group with alfred's summary statistics
(mapping, error rate, insert size, coverage of the kit's targets) and the fraction of targets
covered at least :data:`COVERAGE_LEVELS` times; the charts show the target coverage, on-target
rate, coverage histogram, insert size, base quality by position and GC content; the chromosome
table lists the mapped reads per contig.
"""

from __future__ import annotations

import gzip
import json
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    Chart,
    Column,
    Summary,
    Table,
    bars,
    concat,
    fraction,
    from_rows,
    ratio,
    thin,
)

#: Coverage levels whose fraction of covered targets the summary table lists
COVERAGE_LEVELS = (10, 20, 30, 50, 100)

#: The columns shown in the summary table, in order; the other statistics follow in the details
SUMMARY_COLUMNS = {
    "MappedFraction": fraction(),
    "MappedProperFraction": fraction(),
    "DuplicateFraction": fraction(),
    "OnTargetFraction": fraction(),
    **{f"TargetsAbove{level}x": fraction() for level in COVERAGE_LEVELS},
    "MedianCoverage": bars(precision=1),
    "ErrorRate": Column(precision=4),
    "MedianReadLength": Column(),
    "MedianInsertSize": Column(),
    "MedianMAPQ": Column(),
}

#: Charts of alfred's metrics: metric id -> (name, title, x column, y column, log y)
CHARTS = {
    "targetCoverage": (
        "target_coverage",
        "Targets above coverage",
        "coverage",
        "fraction of targets",
        False,
    ),
    "onTarget": (
        "on_target",
        "On-target rate by target extension",
        "target extension (bp)",
        "fraction on target",
        False,
    ),
    "coverageHistogram": ("coverage", "Coverage distribution", "coverage", "bases", True),
    "insertSize": ("insert_size", "Insert size distribution", "insert size", "read pairs", False),
    "baseQuality": (
        "base_quality",
        "Mean base quality by position",
        "position in read",
        "mean base quality",
        False,
    ),
    "gcContent": ("gc_content", "GC content", "GC fraction", "normalized fraction", False),
}

#: Metrics whose distributions are thinned before plotting
THINNED = ("coverageHistogram", "insertSize")

#: Columns of the chromosome table renamed for the TSV, whose header must not start with ``#``
CHROMOSOME_COLUMNS: dict[str, str] = {"#N": "N", "#GC": "GC"}


def load(path: str) -> dict:
    with gzip.open(path, "rt") as f:
        return json.load(f)


def series(metric: dict, x: str, y: str) -> dict[str | None, pl.DataFrame]:
    """The data series of a metric, by series title, as frames with columns ``x`` and ``y``."""
    xs = metric["x"]["data"][0]["values"]
    return {
        entry.get("title"): pl.DataFrame({x: xs[: len(entry["values"])], y: entry["values"]})
        for entry in metric["y"]["data"]
    }


def _summary_row(label: str, columns: list[str], values: list[Any]) -> dict[str, Any]:
    row: dict[str, Any] = {"library": label}
    for column, value in zip(columns, values):
        if column == "Sample":
            continue
        row["read_group" if column == "Library" else column.lstrip("#")] = value
    row["OnTargetFraction"] = ratio(row.get("AlignedBasesInBed"), row.get("AlignedBases"))
    return row


def _targets_above(metric: dict) -> dict[str, float]:
    frame = next(iter(series(metric, "x", "y").values()))
    by_level = dict(zip(frame["x"].to_list(), frame["y"].to_list()))
    return {
        f"TargetsAbove{level}x": by_level[level] for level in COVERAGE_LEVELS if level in by_level
    }


def _chart_frame(metric_id: str, metric: dict, name: str) -> pl.DataFrame:
    """The rows of one library for a chart: the metric's series in long form."""
    _, _, x, y, _ = CHARTS[metric_id]
    frames = series(metric, x, y)
    library = pl.lit(name).alias("library")
    if metric_id == "insertSize":  # one series per orientation; the sum is the distribution
        total = pl.concat(frames.values()).group_by(x, maintain_order=True).agg(pl.col(y).sum())
        return total.with_columns(library)
    if metric_id == "baseQuality":
        return pl.concat(
            frame.with_columns(library, pl.lit(title).alias("read"))
            for title, frame in frames.items()
        )
    if metric_id == "gcContent":
        return frames["Sample"].with_columns(library)
    return next(iter(frames.values())).with_columns(library)


def summarize(files: dict[str, list[str]], labels: list[str]) -> Summary:
    rows: list[dict[str, Any]] = []
    charts: dict[str, list[pl.DataFrame]] = {metric_id: [] for metric_id in CHARTS}
    chromosomes: list[pl.DataFrame] = []
    reference_curves: dict[str, pl.DataFrame] = {}
    for label, path in zip(labels, files["json"]):
        for sample in load(path)["samples"]:
            summary = sample["summary"]["data"]
            metrics = {
                (group["id"], metric["id"]): metric
                for group in sample["readGroups"]
                for metric in group["metrics"]
            }
            several = len(sample["readGroups"]) > 1
            for values in summary["rows"]:
                row = _summary_row(label, summary["columns"], values)
                group = row["read_group"]
                name = f"{label} ({group})" if several else label
                if target := metrics.get((group, "targetCoverage")):
                    row |= _targets_above(target)
                rows.append(row)
                for metric_id in CHARTS:
                    if metric := metrics.get((group, metric_id)):
                        charts[metric_id].append(_chart_frame(metric_id, metric, name))
                if metric := metrics.get((group, "gcContent")):
                    _, _, x, y, _ = CHARTS["gcContent"]
                    for title, frame in series(metric, x, y).items():
                        if title in ("Reference", "Target") and title not in reference_curves:
                            reference_curves[title] = frame.with_columns(
                                pl.lit(title.lower()).alias("library")
                            )
                if table := metrics.get((group, "mappingByChromosome")):
                    columns = [CHROMOSOME_COLUMNS.get(c, c) for c in table["data"]["columns"]]
                    frame = pl.DataFrame(table["data"]["rows"], schema=columns, orient="row")
                    chromosomes.append(frame.with_columns(pl.lit(name).alias("library")))
    charts["gcContent"] += list(reference_curves.values())

    summary_frame = from_rows(rows, ["library", "read_group", *SUMMARY_COLUMNS])
    rest = [
        c for c in summary_frame.columns if c not in ("library", "read_group", *SUMMARY_COLUMNS)
    ]
    columns = SUMMARY_COLUMNS | {column: Column(detail=True) for column in rest}
    columns["read_group"] = Column(detail=True)

    views: list[Table | Chart] = []
    for metric_id, (name, title, x, y, log_y) in CHARTS.items():
        frame_columns = (
            ["library", "read", x, y] if metric_id == "baseQuality" else ["library", x, y]
        )
        frame = concat(charts[metric_id], frame_columns)
        if metric_id in THINNED:
            frame = thin(frame, x, y)
        views.append(
            Chart(
                name,
                title,
                frame,
                x=x,
                y=y,
                facet="read" if metric_id == "baseQuality" else None,
                log_y=log_y,
            )
        )
    views.append(
        Table(
            "chromosomes",
            "Reads per chromosome",
            concat(chromosomes, ["library", "Chr", "Size", "N", "GC", "GC-fraction", "mapped"]),
            desc="Mapped reads per chromosome, with its size, N and GC content.",
            columns={"mapped": bars()},
            pinned=("library", "Chr"),
        )
    )
    return Summary(
        table=Table(
            "summary",
            "Target coverage",
            summary_frame,
            desc="`alfred qc` of each library: fractions of mapped, properly paired, duplicated "
            "and on-target reads, the fraction of the kit's targets covered at least 10 to 100 "
            "times, median coverage, error rate, read length, insert size and MAPQ. The details "
            "of a row hold the other statistics.",
            columns=columns,
        ),
        views=views,
    )
