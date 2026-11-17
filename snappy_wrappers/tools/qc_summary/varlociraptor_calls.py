# -*- coding: utf-8 -*-
"""The calls of varlociraptor in one VCF, as a table for the report.

One row per call with the posterior probability of every event of the scenario and, per sample,
the allele frequency, the depth, the genotype and the observations. The observations are drawn
in the cell: the reads for the alternative allele to the right and those for the reference (or
a third allele) to the left, one bar per strength of evidence (``SAOBS`` and ``SROBS``), and the
posterior density of the allele frequency (``AFD``) as a bar per sampled frequency. The biases that
varlociraptor found (strand, read orientation, read position, softclip, homopolymer, alt locus,
fragment end) are listed as text. The full observation strings (``OBS``) stay in the details.
"""

from __future__ import annotations

import gzip
import json
import math
import re
from dataclasses import dataclass, field
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    INK,
    TEXT,
    Column,
    Summary,
    Table,
    bars,
    ticks,
)

#: The strengths of evidence of an observation (extended Kass-Raftery scores), weakest first
STRENGTHS = {"E": "equal", "B": "barely", "P": "positive", "S": "strong", "V": "very strong"}

#: The bias estimates of a sample: FORMAT field -> (label, the symbols that mean a bias)
BIASES = {
    "SB": ("strand", "+-"),
    "ROB": ("read orientation", "><"),
    "RPB": ("read position", "^"),
    "SCB": ("softclip", "$"),
    "HE": ("homopolymer", "*"),
    "ALB": ("alt locus", "*"),
    "FEB": ("fragment end", "F"),
}

#: Events whose probability counts against the call, drawn in red
AGAINST = ("ARTIFACT", "ABSENT")

#: The in-cell plot of the observations: bars of the counts per strength, alt to the right
OBSERVATIONS_PLOT = {
    "data": "function(value, row) { return value ? JSON.parse(value) : []; }",
    "spec": {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "width": 110,
        "height": 48,
        "mark": {"type": "bar", "tooltip": True},
        "encoding": {
            "y": {
                "field": "strength",
                "type": "ordinal",
                "sort": list(reversed(STRENGTHS.values())),  # strongest evidence on top
                "axis": {"title": None, "labelFontSize": 8, "ticks": False, "domain": False},
            },
            "x": {
                "field": "reads",
                "type": "quantitative",
                "axis": {"title": None, "labelFontSize": 8, "tickCount": 3, "format": "~s"},
            },
            "color": {
                "field": "allele",
                "type": "nominal",
                "scale": {"domain": ["alt", "ref"], "range": [INK, "#c3c2b7"]},
                "legend": None,
            },
            "opacity": {
                "field": "mapping",
                "type": "nominal",
                "scale": {"domain": ["certain", "uncertain"], "range": [1, 0.4]},
                "legend": None,
            },
            "tooltip": [
                {"field": "allele", "type": "nominal"},
                {"field": "strength", "type": "nominal"},
                {"field": "count", "type": "quantitative"},
                {"field": "mapping", "type": "nominal"},
            ],
        },
        "config": {"view": {"stroke": None}, "axis": {"grid": False, "labelColor": TEXT}},
    },
}

#: The in-cell plot of the allele frequency density: an area over the allele frequency
AFD_PLOT = {
    "data": "function(value, row) { return value ? JSON.parse(value) : []; }",
    "spec": {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "width": 110,
        "height": 36,
        # One bar per sampled allele frequency: three in the discrete case (0, 0.5, 1), a
        # dense row in the continuous one, which reads as the density curve
        "mark": {"type": "bar", "tooltip": True, "color": INK, "opacity": 0.8, "size": 3},
        "encoding": {
            "x": {
                "field": "af",
                "type": "quantitative",
                "scale": {"domain": [0, 1]},
                "axis": {"title": None, "labelFontSize": 8, "values": [0, 0.5, 1]},
            },
            "y": {"field": "density", "type": "quantitative", "axis": None},
            "tooltip": [
                {"field": "af", "type": "quantitative", "title": "allele frequency"},
                {"field": "density", "type": "quantitative", "format": ".3~g"},
            ],
        },
        "config": {"view": {"stroke": None}, "axis": {"grid": False, "labelColor": TEXT}},
    },
}


@dataclass
class Calls:
    """The records of a varlociraptor VCF: its samples, events and one dict per call."""

    samples: list[str]
    events: list[str]
    rows: list[dict[str, Any]] = field(default_factory=list)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def phred_to_probability(text: str) -> float | None:
    """A probability from varlociraptor's PHRED-scaled value; ``inf`` is zero, ``.`` unknown."""
    if text in (".", ""):
        return None
    value = float(text)
    if math.isinf(value):
        return 0.0
    return 10 ** (-value / 10)


def observations(saobs: str, srobs: str) -> list[dict[str, Any]]:
    """The simplified observations of a sample as records: one per allele (``alt`` for the
    reads favoring the alternative allele, ``ref`` for those favoring the reference or a third
    allele), strength and mapping certainty, with the count and the signed ``reads`` of the
    plot (alt to the right, ref to the left)."""
    records = []
    for allele, text, sign in (("alt", saobs, 1), ("ref", srobs, -1)):
        for count, code in re.findall(r"(\d+)([A-Za-z])", text or ""):
            strength = STRENGTHS.get(code.upper())
            if strength is None:
                continue
            records.append(
                {
                    "allele": allele,
                    "strength": strength,
                    "mapping": "certain" if code.isupper() else "uncertain",
                    "count": int(count),
                    "reads": sign * int(count),
                }
            )
    return records


def allele_frequency_density(afd: str) -> list[dict[str, float]]:
    """The posterior density of the allele frequency from ``AFD`` (``af=phred,...``)."""
    records = []
    for entry in (afd or "").split(","):
        if "=" not in entry:
            continue
        af, phred = entry.split("=", 1)
        records.append({"af": float(af), "density": 10 ** (-float(phred) / 10)})
    return records


def biases(sample: dict[str, str]) -> str:
    """The biases that varlociraptor found for a sample, as text; empty without any."""
    found = []
    for key, (label, symbols) in BIASES.items():
        value = sample.get(key, ".")
        if value and any(symbol in value for symbol in symbols):
            found.append(f"{label} {value}")
    return ", ".join(found)


def variant_type(ref: str, alt: str, info: dict[str, str]) -> str:
    if "SVTYPE" in info:
        return info["SVTYPE"]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INS" if len(alt) > len(ref) else "DEL"


def read_calls(path: str) -> Calls:
    """Read the calls of a varlociraptor VCF (one ALT per record, as varlociraptor writes)."""
    samples: list[str] = []
    events: list[str] = []
    rows: list[dict[str, Any]] = []
    with _open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                if match := re.match(r"##INFO=<ID=PROB_(\w+),", line):
                    events.append(match.group(1))
                continue
            if line.startswith("#"):
                samples = line.split("\t")[9:]
                continue
            fields = line.split("\t")
            chromosome, position, _, ref, alt = fields[:5]
            info: dict[str, str] = {}
            for entry in fields[7].split(";"):
                key, _, value = entry.partition("=")
                info[key] = value
            row: dict[str, Any] = {
                "chromosome": chromosome,
                "position": int(position),
                "type": variant_type(ref, alt, info),
                "ref": ref,
                "alt": alt,
            }
            for event in events:
                row[f"prob: {event.lower()}"] = phred_to_probability(info.get(f"PROB_{event}", "."))
            keys = fields[8].split(":")
            for sample, text in zip(samples, fields[9:]):
                values = dict(zip(keys, text.split(":")))
                af = values.get("AF", ".")
                dp = values.get("DP", ".")
                row[f"{sample}: allele frequency"] = None if af == "." else float(af)
                row[f"{sample}: depth"] = None if dp == "." else int(dp)
                row[f"{sample}: genotype"] = values.get("GT", ".")
                row[f"{sample}: observations"] = json.dumps(
                    observations(values.get("SAOBS", ""), values.get("SROBS", ""))
                )
                row[f"{sample}: allele frequency density"] = json.dumps(
                    allele_frequency_density(values.get("AFD", ""))
                )
                row[f"{sample}: biases"] = biases(values)
                row[f"{sample}: raw observations"] = values.get("OBS", "")
            row["hints"] = info.get("HINTS", "")
            rows.append(row)
    return Calls(samples, events, rows)


def _probability(event: str) -> Column:
    scheme = "reds" if event.upper() in AGAINST else "blues"
    plot = {
        "heatmap": {"scale": "linear", "domain": ["0", "1"], "color-scheme": scheme},
    }
    return Column(precision=3, plot=plot)


def summarize(path: str) -> Summary:
    """The calls table of one VCF, ready to render."""
    calls = read_calls(path)
    columns: dict[str, Column] = {"position": Column(precision=0), "hints": Column(detail=True)}
    for event in calls.events:
        columns[f"prob: {event.lower()}"] = _probability(event)
    for sample in calls.samples:
        columns[f"{sample}: allele frequency"] = ticks(precision=3)
        columns[f"{sample}: depth"] = bars()
        columns[f"{sample}: genotype"] = Column()
        columns[f"{sample}: observations"] = Column(custom_plot=OBSERVATIONS_PLOT)
        columns[f"{sample}: allele frequency density"] = Column(custom_plot=AFD_PLOT)
        columns[f"{sample}: biases"] = Column()
        columns[f"{sample}: raw observations"] = Column(detail=True)
    order = ["chromosome", "position", "type", "ref", "alt"]
    frame = (
        pl.from_dicts(calls.rows, infer_schema_length=None)
        if calls.rows
        else pl.DataFrame(schema={column: pl.Null for column in [*order, *columns]})
    )
    rest = [c for c in frame.columns if c not in order]
    frame = frame.select(*order, *rest)
    return Summary(
        table=Table(
            "calls",
            "varlociraptor calls",
            frame,
            desc="One row per call. `prob:` columns are the posterior probabilities of the "
            "events of the scenario (red for artifact and absent). Per sample: the allele "
            "frequency with its posterior density, the depth, the genotype, the observations "
            "(bars per strength of evidence: reads for the alternative allele to the right, "
            "for the reference or another allele to the left; pale bars from reads with "
            "uncertain mapping) and the biases varlociraptor found. The details hold the raw "
            "observation strings and the hints.",
            columns=columns,
            pinned=("chromosome", "position"),
            page_size=20,
        ),
        views=[],
    )
