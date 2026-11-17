# -*- coding: utf-8 -*-
"""The ``bcftools stats`` files of a ``variant_calling`` task, one per sample.

The summary table has one row per sample and variant set (``known`` and ``novel`` with
``--split-by-ID``, else ``all``) with the per-sample counts (``PSC``), the Ts/Tv ratio and the
site counts (``SN``); the charts show the substitution types, the indel length, depth and quality
distributions (``ST``, ``IDD``, ``DP``, ``QUAL``).
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
    from_rows,
    numeric,
    ratio,
    thin,
)

#: ``PSC`` columns of the summary table and their names
PSC_COLUMNS = {
    "nRefHom": "hom_ref",
    "nNonRefHom": "hom_alt",
    "nHets": "het",
    "nTransitions": "transitions",
    "nTransversions": "transversions",
    "nIndels": "indels",
    "average depth": "average_depth",
    "nSingletons": "singletons",
    "nMissing": "missing",
}

#: ``SN`` keys of the summary table and their names
SN_COLUMNS = {
    "number of records": "sites",
    "number of SNPs": "snv_sites",
    "number of MNPs": "mnv_sites",
    "number of indels": "indel_sites",
    "number of others": "other_sites",
    "number of multiallelic sites": "multiallelic_sites",
}

SUMMARY_COLUMNS = {
    "het": bars(),
    "hom_alt": bars(),
    "het_hom_ratio": Column(precision=2),
    "ts_tv": Column(precision=2),
    "indels": bars(),
    "singletons": bars(),
    "missing": bars(),
    "average_depth": bars(precision=1),
    "sites": Column(detail=True),
    "snv_sites": Column(detail=True),
    "mnv_sites": Column(detail=True),
    "indel_sites": Column(detail=True),
    "other_sites": Column(detail=True),
    "multiallelic_sites": Column(detail=True),
    "hom_ref": Column(detail=True),
    "transitions": Column(detail=True),
    "transversions": Column(detail=True),
}

#: Charts of the distribution sections: section -> (name, title, x column, x title, y column)
CHARTS = {
    "ST": ("substitutions", "Substitution types", "type", "substitution", "count"),
    "IDD": (
        "indel_length",
        "Indel length distribution",
        "length (deletions negative)",
        "length (deletions negative)",
        "number of sites",
    ),
    "DP": ("depth", "Depth distribution", "bin", "depth", "fraction of genotypes (%)"),
    "QUAL": ("quality", "Quality distribution", "Quality", "quality", "number of SNPs"),
}


def parse(path: str) -> tuple[dict[str, pl.DataFrame], dict[str, str]]:
    """Return the sections of a ``bcftools stats`` file as frames (columns named as in the
    header, without the ``[n]`` prefixes) and the label of each variant set id."""
    headers: dict[str, list[str]] = {}
    rows: dict[str, list[list[str]]] = {}
    labels: dict[str, str] = {}
    with open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if match := re.match(r"# (\w+)\t\[2\]id\t(.*)", line):
                section, columns = match.groups()
                headers[section] = [
                    "id",
                    *(re.sub(r"^\[\d+\]", "", c) for c in columns.split("\t")),
                ]
            elif line.startswith("ID\t"):
                _, set_id, text = line.split("\t", 2)
                labels[set_id] = _set_label(text)
            elif line and not line.startswith("#"):
                section, *fields = line.split("\t")
                rows.setdefault(section, []).append(fields)
    frames = {}
    for section, header in headers.items():
        frame = pl.DataFrame(
            [row[: len(header)] for row in rows.get(section, [])], schema=header, orient="row"
        )
        # Keys, sample names and ``>500`` depth bins stay text
        frames[section] = numeric(frame, [c for c in header if c != "id"])
    return frames, labels


def _set_label(text: str) -> str:
    """``known`` or ``novel`` after the file name of a ``--split-by-ID`` run, else ``all``."""
    if match := re.search(r":(\w+) \(", text):
        return match.group(1)
    return "all"


def _with_set(frame: pl.DataFrame, sets: dict[str, str], label: str) -> pl.DataFrame:
    return frame.with_columns(
        pl.lit(label).alias("library"),
        pl.col("id").replace_strict(sets, default="all").alias("set"),
    )


def summarize(files: dict[str, list[str]], labels: list[str]) -> Summary:
    rows: list[dict[str, Any]] = []
    charts: dict[str, list[pl.DataFrame]] = {section: [] for section in CHARTS}
    for label, path in zip(labels, files["txt"]):
        frames, sets = parse(path)
        numbers = {
            (row["id"], row["key"].rstrip(":")): row["value"]
            for row in frames.get("SN", pl.DataFrame()).iter_rows(named=True)
        }
        for psc in frames.get("PSC", pl.DataFrame()).iter_rows(named=True):
            # The task writes one file per sample; a file of several samples has a row each
            row = {"library": psc.get("sample", label), "set": sets.get(psc["id"], "all")}
            row |= {name: psc.get(key) for key, name in PSC_COLUMNS.items()}
            row |= {name: numbers.get((psc["id"], key)) for key, name in SN_COLUMNS.items()}
            row["ts_tv"] = ratio(row.get("transitions"), row.get("transversions"))
            row["het_hom_ratio"] = ratio(row.get("het"), row.get("hom_alt"))
            rows.append(row)
        for section, (_, _, x, _, y) in CHARTS.items():
            if (frame := frames.get(section)) is not None and not frame.is_empty():
                charts[section].append(_with_set(frame, sets, label).select("library", "set", x, y))

    summary = from_rows(rows, ["library", "set", *SUMMARY_COLUMNS])
    facet = "set" if summary.height and summary["set"].n_unique() > 1 else None
    views: list[Table | Chart] = []
    for section, (name, title, x, x_title, y) in CHARTS.items():
        frame = concat(charts[section], ["library", "set", x, y])
        if section == "DP":  # the last bin is ``>500``, which has no place on the axis
            frame = frame.with_columns(pl.col(x).cast(pl.Float64, strict=False)).drop_nulls(x)
        if section == "QUAL":
            frame = thin(frame, x, y, by=("library", "set"))
        views.append(
            Chart(
                name,
                title,
                frame,
                x=x,
                y=y,
                x_title=x_title,
                x_type="nominal" if section == "ST" else "quantitative",
                mark="bar" if section == "ST" else "line",
                facet=facet,
            )
        )
    return Summary(
        table=Table(
            "summary",
            "Variant statistics",
            summary,
            desc="`bcftools stats` of each sample: genotype counts, het/hom and Ts/Tv ratios, "
            "indels, singletons, missing genotypes and depth; `known` and `novel` sets by the "
            "ID column. The details of a row hold the site counts of the VCF.",
            columns=SUMMARY_COLUMNS,
            pinned=("library", "set"),
        ),
        views=views,
    )
