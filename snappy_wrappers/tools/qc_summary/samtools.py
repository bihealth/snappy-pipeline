# -*- coding: utf-8 -*-
"""The ``samtools stats``, ``flagstat`` and ``idxstats`` files of a mapping task.

The summary table has one row per library with the counts and fractions of the ``SN`` section
and the singletons of ``flagstat``; the charts show the insert size and coverage distributions
and the mean base quality by cycle (``IS``, ``COV``, ``FFQ``, ``LFQ``); the idxstats table lists
the mapped reads per contig.
"""

from __future__ import annotations

import re
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
    number,
    ratio,
    thin,
)

#: ``SN`` keys of the summary table and their column names
SN_COLUMNS = {
    "raw total sequences": "reads",
    "reads mapped": "mapped",
    "reads mapped and paired": "mapped_and_paired",
    "reads properly paired": "properly_paired",
    "reads duplicated": "duplicated",
    "reads MQ0": "mapq0",
    "supplementary alignments": "supplementary",
    "non-primary alignments": "secondary",
    "pairs on different chromosomes": "pairs_on_different_chromosomes",
    "error rate": "error_rate",
    "average length": "average_length",
    "average quality": "average_quality",
    "insert size average": "insert_size_average",
    "insert size standard deviation": "insert_size_sd",
    "bases mapped (cigar)": "bases_mapped",
}

#: ``flagstat`` lines of the summary table and their column names
FLAGSTAT_COLUMNS = {
    "singletons": "singletons",
    "with mate mapped to a different chr (mapQ>=5)": "mate_on_other_chromosome",
}

#: The columns shown in the summary table, in order; the counts follow in the row details
SUMMARY_COLUMNS = {
    "reads": bars(),
    "mapped_fraction": fraction(),
    "properly_paired_fraction": fraction(),
    "duplicated_fraction": fraction(),
    "mapq0_fraction": fraction(),
    "error_rate": Column(precision=4),
    "average_length": Column(precision=0),
    "average_quality": Column(precision=1),
    "insert_size_average": Column(precision=0),
    "insert_size_sd": Column(precision=0),
}

COLUMNS = SUMMARY_COLUMNS | {
    column: Column(detail=True)
    for column in (*SN_COLUMNS.values(), *FLAGSTAT_COLUMNS.values())
    if column not in SUMMARY_COLUMNS
}


def parse_stats(path: str) -> dict[str, list[list[str]]]:
    """Return the rows of each section of a ``samtools stats`` file, by section name."""
    sections: dict[str, list[list[str]]] = {}
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            sections.setdefault(fields[0], []).append(fields[1:])
    return sections


def summary_numbers(sections: dict[str, list[list[str]]]) -> dict[str, float | int]:
    """The ``SN`` section as key -> number (the key without its colon)."""
    return {row[0].rstrip(":"): number(row[1]) for row in sections.get("SN", ())}


def parse_flagstat(path: str) -> dict[str, int]:
    """The QC-passed count of each ``flagstat`` line, by its label without the percentages."""
    result = {}
    with open(path, "rt") as f:
        for line in f:
            if match := re.match(r"(\d+) \+ (\d+) (.+)", line.strip()):
                label = re.sub(r" \(.*\)$", "", match.group(3))
                result[label] = int(match.group(1))
    return result


def parse_idxstats(path: str) -> pl.DataFrame:
    """The ``idxstats`` table of the placed contigs: contig, length, mapped, unmapped."""
    frame = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["contig", "length", "mapped", "unmapped"],
        schema_overrides={"contig": pl.String},
    )
    return frame.filter(pl.col("contig") != "*")


def _distribution(sections, section: str, x: str, y: str, columns: tuple[int, int]):
    rows = sections.get(section, ())
    return pl.DataFrame(
        {
            x: [number(row[columns[0]]) for row in rows],
            y: [number(row[columns[1]]) for row in rows],
        }
    )


def _quality_by_cycle(sections, section: str) -> pl.DataFrame:
    """Mean base quality per cycle from the counts per quality of ``FFQ`` or ``LFQ``."""
    cycles, means = [], []
    for row in sections.get(section, ()):
        counts = [int(value) for value in row[1:]]
        total = sum(counts)
        if total:
            cycles.append(int(row[0]))
            means.append(sum(q * n for q, n in enumerate(counts)) / total)
    return pl.DataFrame({"cycle": cycles, "mean_quality": means})


def summarize(files: dict[str, list[str]], labels: list[str]) -> Summary:
    rows, insert_sizes, coverages, qualities, idxstats = [], [], [], [], []
    for label, stats_path, flagstat_path, idxstats_path in zip(
        labels, files["bamstats"], files["flagstats"], files["idxstats"]
    ):
        sections = parse_stats(stats_path)
        numbers = summary_numbers(sections)
        flagstat = parse_flagstat(flagstat_path)
        row: dict[str, Any] = {"library": label}
        row |= {column: numbers.get(key) for key, column in SN_COLUMNS.items()}
        row |= {column: flagstat.get(key) for key, column in FLAGSTAT_COLUMNS.items()}
        reads, mapped = numbers.get("raw total sequences"), numbers.get("reads mapped")
        row["mapped_fraction"] = ratio(mapped, reads)
        row["properly_paired_fraction"] = ratio(
            numbers.get("reads properly paired"), numbers.get("reads paired")
        )
        row["duplicated_fraction"] = ratio(numbers.get("reads duplicated"), reads)
        row["mapq0_fraction"] = ratio(numbers.get("reads MQ0"), mapped)
        rows.append(row)
        library = pl.lit(label).alias("library")
        insert_sizes.append(
            _distribution(sections, "IS", "insert_size", "pairs", (0, 1)).with_columns(library)
        )
        coverages.append(
            _distribution(sections, "COV", "coverage", "bases", (1, 2)).with_columns(library)
        )
        for read, section in (("first", "FFQ"), ("last", "LFQ")):
            qualities.append(
                _quality_by_cycle(sections, section).with_columns(
                    library, pl.lit(read).alias("read")
                )
            )
        idxstats.append(parse_idxstats(idxstats_path).with_columns(library))

    summary = from_rows(rows, ["library", *COLUMNS])
    contigs = concat(idxstats, ["library", "contig", "length", "mapped", "unmapped"])
    contigs = contigs.with_columns(
        (pl.col("mapped") / pl.col("length") * 1e6).alias("mapped_per_mb")
    )
    mapped_per_mb = contigs["mapped_per_mb"].max() if contigs.height else None
    return Summary(
        table=Table(
            "summary",
            "BAM QC",
            summary,
            desc="`samtools stats` of each library: read counts, fractions of mapped, properly "
            "paired, duplicated and MAPQ 0 reads, error rate, read length, base quality and "
            "insert size. The details of a row hold the counts.",
            columns=COLUMNS,
        ),
        views=[
            Chart(
                "insert_size",
                "Insert size distribution",
                thin(
                    concat(insert_sizes, ["library", "insert_size", "pairs"]),
                    "insert_size",
                    "pairs",
                ),
                x="insert_size",
                y="pairs",
                x_title="insert size",
                y_title="read pairs",
            ),
            Chart(
                "coverage",
                "Coverage distribution",
                thin(concat(coverages, ["library", "coverage", "bases"]), "coverage", "bases"),
                x="coverage",
                y="bases",
                y_title="bases of the reference",
                log_y=True,
            ),
            Chart(
                "base_quality",
                "Mean base quality by cycle",
                concat(qualities, ["library", "read", "cycle", "mean_quality"]),
                x="cycle",
                y="mean_quality",
                y_title="mean base quality",
                facet="read",
            ),
            Table(
                "idxstats",
                "Reads per contig",
                contigs,
                desc="`samtools idxstats`: the mapped and unmapped reads placed on each contig.",
                columns={
                    "mapped": bars(),
                    "mapped_per_mb": fraction(
                        label="mapped per Mb", high=float(mapped_per_mb or 1.0)
                    ),
                },
                pinned=("library", "contig"),
            ),
        ],
    )
