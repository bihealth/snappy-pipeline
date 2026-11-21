# -*- coding: utf-8 -*-
"""The columns of varlociraptor's calls:
the posterior probability of every event of the scenario and, per sample, the allele frequency,
the depth, the genotype and the observations.

The observations are drawn in the cell:
the reads for the alternative allele to the right and those for the reference (or a third allele)
to the left, one bar per strength of evidence (``SAOBS`` and ``SROBS``),
and the posterior density of the allele frequency (``AFD``) as a bar per sampled frequency,
relative to the cell's maximum.
The cells hold the compact strings of the VCF (``SAOBS|SROBS``,
and ``AFD`` thinned to :data:`AFD_POINTS` frequencies), which the browser parses,
so a row costs tens of bytes per sample.
The biases that varlociraptor found (strand, read orientation, read position, softclip, homopolymer,
alt locus, fragment end) are listed as text.
"""

from __future__ import annotations

import math
import re
from typing import Any
from urllib.parse import quote, urlencode

from snappy_wrappers.tools.qc_summary import (
    AGAINST,
    CELL_CONFIG,
    INK,
    MUTED,
    SEQUENTIAL,
    VEGA_LITE,
    Column,
    bars,
    heatmap,
    ticks,
)

from .vcf import Record, Vcf, short_labels

#: What the table's description says about the varlociraptor columns
COLUMN_NOTES = [
    "`event`: the most probable event of the call; `P(...)`: the posterior probability of "
    "every event of the scenario, artifact and absent on the red ramp.",
    "Per sample: `AF` (allele frequency), `DP` (depth) and `GT` (genotype).",
    "Row details, per sample: the observations as bars per strength of evidence, reads for the "
    "alternative allele to the right, for the reference or another allele to the left, pale "
    "from reads with uncertain mapping; and the posterior density of the allele frequency, "
    "relative to the cell's maximum.",
    "`biases` (on demand): the biases varlociraptor found.",
]

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

#: varlociraptor-inspect, which draws the observations of a record from the URL's parameters
INSPECT_URL = "https://varlociraptor.github.io/varlociraptor-inspect/"

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

#: The JS function that turns a cell ``af=phred,...`` into the records of the density plot,
#: the density relative to the cell's maximum
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

    Every strength has its row, filled or not,
    and the x axis spans ``max_reads`` both ways (the most reads of one strength in the column),
    so the cells of a column share their scales and the rows compare.
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


def phred_to_probability(text: str) -> float | None:
    """A probability from varlociraptor's PHRED-scaled value; ``inf`` is zero, ``.`` unknown."""
    if text in (".", ""):
        return None
    value = float(text)
    if math.isinf(value):
        return 0.0
    return 10 ** (-value / 10)


def observations(saobs: str, srobs: str) -> list[dict[str, Any]]:
    """The simplified observations of a sample as records:
    one per allele (``alt`` for the reads favoring the alternative allele,
    ``ref`` for those favoring the reference or a third allele), strength and mapping certainty,
    with the count and the signed ``reads`` of the plot (alt to the right, ref to the left)."""
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


def events(vcf: Vcf) -> list[str]:
    """The events of the scenario, from the ``PROB_*`` INFO fields of the header."""
    return [name[5:] for name in vcf.info if name.startswith("PROB_")]


def probability_column(event: str) -> str:
    return f"P({event.lower()})"


def probabilities(record: Record, events: list[str]) -> dict[str, Any]:
    """One posterior probability per event."""
    return {
        probability_column(event): phred_to_probability(record.info.get(f"PROB_{event}", "."))
        for event in events
    }


def sample_values(record: Record, samples: list[str]) -> dict[str, Any]:
    """Per sample the allele frequency, depth, genotype,
    the compact observation and density strings and the biases."""
    row: dict[str, Any] = {}
    for sample in samples:
        fields = record.samples.get(sample, {})
        af = fields.get("AF", ".")
        dp = fields.get("DP", ".")
        row[f"{sample}: AF"] = None if af == "." else float(af)
        row[f"{sample}: DP"] = None if dp == "." else int(dp)
        row[f"{sample}: GT"] = fields.get("GT", ".")
        saobs, srobs = fields.get("SAOBS", ""), fields.get("SROBS", "")
        row[f"{sample}: observations"] = f"{saobs}|{srobs}" if saobs or srobs else ""
        row[f"{sample}: AFD"] = compact_afd(fields.get("AFD", ""))
        row[f"{sample}: biases"] = biases(fields)
    return row


def inspect_url(record: Record, events: list[str], samples: list[str]) -> str:
    """The link that opens the record in varlociraptor-inspect:
    the event probabilities and per sample the observations and the allele frequency density, nothing else.
    The samples are ``s1``, ``s2``, ...
    in the order of the VCF, so the URL carries no names."""
    params = [(f"PROB_{event}", record.info.get(f"PROB_{event}", "")) for event in events]
    for number, sample in enumerate(samples, 1):
        fields = record.samples.get(sample, {})
        params.append((f"OBS_s{number}", fields.get("OBS", "")))
        params.append((f"AFD_s{number}", fields.get("AFD", "")))
    return INSPECT_URL + "?" + urlencode(params, quote_via=quote)


def inspect_values(record: Record, events: list[str], samples: list[str]) -> dict[str, Any]:
    return {"inspect": "inspect", "inspect_url": inspect_url(record, events, samples)}


def inspect_columns() -> dict[str, Column]:
    return {
        "inspect": Column(links={"varlociraptor-inspect": "{inspect_url}"}),
        "inspect_url": Column(display="hidden"),
    }


def probability_columns(vcf: Vcf) -> dict[str, Column]:
    """The probabilities as heatmaps, the against-events on the red ramp."""
    result: dict[str, Column] = {}
    for event in events(vcf):
        colors = AGAINST if event.upper() in AGAINST_EVENTS else SEQUENTIAL
        result[probability_column(event)] = heatmap(colors=colors)
    return result


def sample_columns(vcf: Vcf, rows: list[dict[str, Any]]) -> dict[str, Column]:
    """The rendering of the per-sample columns, in order, labelled with the short sample names."""
    result: dict[str, Column] = {}
    labels = short_labels(vcf.samples)
    for sample in vcf.samples:
        short = labels[sample]
        result[f"{sample}: AF"] = ticks(label=f"{short}: AF", precision=3, domain=(0, 1))
        result[f"{sample}: DP"] = bars(label=f"{short}: DP")
        result[f"{sample}: GT"] = Column(label=f"{short}: GT")
        most = max(
            (
                record["count"]
                for row in rows
                if row[f"{sample}: observations"]
                for record in observations(*row[f"{sample}: observations"].split("|"))
            ),
            default=0,
        )
        result[f"{sample}: observations"] = Column(
            label=f"{short}: observations", custom_plot=observations_plot(most), display="detail"
        )
        result[f"{sample}: AFD"] = Column(
            label=f"{short}: AFD", custom_plot=AFD_PLOT, display="detail"
        )
        result[f"{sample}: biases"] = Column(label=f"{short}: biases", display="available")
    return result
