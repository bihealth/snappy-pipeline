# -*- coding: utf-8 -*-
"""The calls of varlociraptor in one VCF, as a table for the report.

One row per call with the posterior probability of every event of the scenario and, per sample,
the allele frequency, the depth, the genotype and the observations.
The observations are drawn in the cell: the reads for the alternative allele to the right and those for the reference (or a third allele) to the left, one bar per strength of evidence (``SAOBS`` and ``SROBS``), and the posterior density of the allele frequency (``AFD``) as a bar per sampled frequency, relative to the cell's maximum.
The cells hold the compact strings of the VCF (``SAOBS|SROBS``, and ``AFD`` thinned to :data:`AFD_POINTS` frequencies), which the browser parses, so a row costs tens of bytes per sample.
The biases that varlociraptor found (strand, read orientation, read position, softclip, homopolymer, alt locus, fragment end) are listed as text.
"""

from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass, field
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    AGAINST,
    CELL_CONFIG,
    INK,
    MUTED,
    SEQUENTIAL,
    VEGA_LITE,
    Column,
    Summary,
    Table,
    bars,
    heatmap,
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

#: Events whose probability counts against the call, drawn on the red ramp
AGAINST_EVENTS = ("ARTIFACT", "ABSENT")

#: Sampled allele frequencies that a density cell keeps at most
AFD_POINTS = 24

#: The JS function that turns a cell ``SAOBS|SROBS`` into the records of the observations plot:
#: one per allele, strength and mapping certainty, the count signed (alt to the right)
OBSERVATIONS_DATA = """function(value, row) {
  const strengths = {E: "equal", B: "barely", P: "positive", S: "strong", V: "very strong"};
  const records = [];
  if (!value) return records;
  const [saobs, srobs] = value.split("|");
  for (const [allele, text, sign] of [["alt", saobs, 1], ["ref", srobs, -1]]) {
    for (const match of (text || "").matchAll(/(\\d+)([A-Za-z])/g)) {
      const strength = strengths[match[2].toUpperCase()];
      if (!strength) continue;
      const count = parseInt(match[1]);
      const mapping = match[2] === match[2].toUpperCase() ? "certain" : "uncertain";
      records.push({allele: allele, strength: strength, mapping: mapping, count: count, reads: sign * count});
    }
  }
  return records;
}"""

#: The JS function that turns a cell ``af=phred,...`` into the records of the density plot, the
#: density relative to the cell's maximum
AFD_DATA = """function(value, row) {
  if (!value) return [];
  const records = value.split(",").filter((entry) => entry.includes("=")).map((entry) => {
    const [af, phred] = entry.split("=");
    return {af: parseFloat(af), density: Math.pow(10, -parseFloat(phred) / 10)};
  });
  const max = Math.max(0, ...records.map((record) => record.density));
  return records.map((record) => ({af: record.af, density: max > 0 ? record.density / max : 0}));
}"""


def observations_plot(max_reads: int) -> dict:
    """The in-cell plot of the observations: bars of the counts per strength, alt to the right.

    Every strength has its row, filled or not, and the x axis spans ``max_reads`` both ways (the
    most reads of one strength in the column), so the cells of a column share their scales and
    the rows compare.
    """
    limit = max(1, max_reads)
    return {
        "data": OBSERVATIONS_DATA,
        "spec": {
            "$schema": VEGA_LITE,
            "width": 110,
            "height": 60,
            "mark": {"type": "bar", "tooltip": True},
            "encoding": {
                "y": {
                    "field": "strength",
                    "type": "ordinal",
                    # Every strength, the strongest evidence on top
                    "scale": {"domain": list(reversed(STRENGTHS.values()))},
                    "axis": {"title": None, "domain": False},
                },
                "x": {
                    "field": "reads",
                    "type": "quantitative",
                    "scale": {"domain": [-limit, limit]},
                    "axis": {"title": None, "values": [-limit, 0, limit], "format": "~s"},
                },
                "color": {
                    "field": "allele",
                    "type": "nominal",
                    "scale": {"domain": ["alt", "ref"], "range": [INK, MUTED]},
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
            "config": CELL_CONFIG,
        },
    }


#: The in-cell plot of the allele frequency density: bars over the allele frequency
AFD_PLOT = {
    "data": AFD_DATA,
    "spec": {
        "$schema": VEGA_LITE,
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
                "axis": {"title": None, "values": [0, 0.5, 1]},
            },
            "y": {"field": "density", "type": "quantitative", "axis": None},
            "tooltip": [
                {"field": "af", "type": "quantitative", "title": "allele frequency"},
                {
                    "field": "density",
                    "type": "quantitative",
                    "format": ".3~g",
                    "title": "relative density",
                },
            ],
        },
        "config": CELL_CONFIG,
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


def compact_afd(afd: str, points: int = AFD_POINTS) -> str:
    """``AFD`` (``af=phred,...``) with at most ``points`` frequencies, the mode always among them,
    and short numbers; entries without density (``inf``) are dropped."""
    entries = []
    for entry in (afd or "").split(","):
        if "=" not in entry:
            continue
        af, phred = entry.split("=", 1)
        if not math.isinf(float(phred)):
            entries.append((float(af), float(phred)))
    if len(entries) > points:
        keep = {round(i * (len(entries) - 1) / (points - 2)) for i in range(points - 1)}
        keep.add(min(range(len(entries)), key=lambda i: entries[i][1]))
        entries = [entries[i] for i in sorted(keep)]
    return ",".join(f"{round(af, 3):g}={round(phred, 1) or 0.0:.1f}" for af, phred in entries)


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
                saobs, srobs = values.get("SAOBS", ""), values.get("SROBS", "")
                row[f"{sample}: observations"] = f"{saobs}|{srobs}" if saobs or srobs else ""
                row[f"{sample}: allele frequency density"] = compact_afd(values.get("AFD", ""))
                row[f"{sample}: biases"] = biases(values)
            row["hints"] = info.get("HINTS", "")
            rows.append(row)
    return Calls(samples, events, rows)


def _probability(event: str) -> Column:
    return heatmap(colors=AGAINST if event.upper() in AGAINST_EVENTS else SEQUENTIAL)


def summarize(path: str) -> Summary:
    """The calls table of one VCF, ready to render."""
    calls = read_calls(path)
    columns: dict[str, Column] = {"position": Column(precision=0), "hints": Column(detail=True)}
    for event in calls.events:
        columns[f"prob: {event.lower()}"] = _probability(event)
    for sample in calls.samples:
        columns[f"{sample}: allele frequency"] = ticks(precision=3, domain=(0, 1))
        columns[f"{sample}: depth"] = bars()
        columns[f"{sample}: genotype"] = Column()
        most = max(
            (
                record["count"]
                for row in calls.rows
                for record in observations(*row[f"{sample}: observations"].split("|"))
                if row[f"{sample}: observations"]
            ),
            default=0,
        )
        columns[f"{sample}: observations"] = Column(custom_plot=observations_plot(most))
        columns[f"{sample}: allele frequency density"] = Column(custom_plot=AFD_PLOT)
        columns[f"{sample}: biases"] = Column()
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
            "frequency with its posterior density (relative to the cell's maximum), the depth, "
            "the genotype, the observations (bars per strength of evidence: reads for the "
            "alternative allele to the right, for the reference or another allele to the left; "
            "pale bars from reads with uncertain mapping) and the biases varlociraptor found. "
            "The details hold the hints.",
            columns=columns,
            pinned=("chromosome", "position"),
            page_size=20,
        ),
        views=[],
    )
