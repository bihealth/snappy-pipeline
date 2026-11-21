# -*- coding: utf-8 -*-
"""The pedigree of a call as a plot in the row details: who carries the variant, and how surely.

The members of the entity's pedigree, from the sample sheet, are laid out once per table:
generations as rows, children below their parents, partners side by side, with the unsequenced donors that connect members and the unknown parents of members with one parent drawn hollow.
Per call, the plot reads the per-sample columns of the row: the allele frequency as a bar under the node with the 90 % credible interval from the allele frequency density, the depth and the genotype as text, and the call's event with its posterior as the title.
The glyphs follow the pedigree convention: square for male, circle for female, diamond for unknown sex, filled when affected.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

from snappy_wrappers.tools.qc_summary import CELL_CONFIG, INK, MUTED, TEXT, VEGA_LITE, Column

#: Pixels per unit of the layout, sideways and downwards
X_STEP, Y_STEP = 96, 104

#: Characters of a member's label under its node; longer labels are cut, the tooltip has them
LABEL_LENGTH = 14

#: The sex of the unknown parent that fills the other slot
OPPOSITE = {"father": "male", "mother": "female"}


@dataclass
class Member:
    """A node of the pedigree: a sequenced member, an unsequenced donor or an unknown parent."""

    name: str
    label: str
    father: str = "0"
    mother: str = "0"
    sex: str = "unknown"
    #: ``True`` affected, ``False`` unaffected, ``None`` unknown
    affected: bool | None = None
    #: The VCF sample of the member, ``None`` without one
    sample: str | None = None


def members_of(pedigree: list[dict[str, Any]], samples: list[str]) -> list[Member]:
    """The members from the params of the report (``name``, ``label``, ``father``, ``mother``,
    ``sex``, ``affected``, ``sample``), the samples that the VCF lacks unobserved, plus an
    unknown parent for every member with one parent among them."""
    members = [
        Member(
            name=entry["name"],
            label=entry.get("label") or entry["name"],
            father=entry.get("father") or "0",
            mother=entry.get("mother") or "0",
            sex=entry.get("sex") or "unknown",
            affected=entry.get("affected"),
            sample=entry.get("sample") if entry.get("sample") in samples else None,
        )
        for entry in pedigree
    ]
    names = {member.name for member in members}
    for member in list(members):
        known = {"father": member.father in names, "mother": member.mother in names}
        if known["father"] == known["mother"]:
            continue
        missing = "mother" if known["father"] else "father"
        sex = OPPOSITE[missing]
        name = f"{member.name}.{missing}"
        setattr(member, missing, name)
        members.append(Member(name=name, label="?", sex=sex))
        names.add(name)
    return members


def generations(members: list[Member]) -> dict[str, int]:
    """Member -> generation from the top: below both parents, and a founder in the row above
    its children (so that partners share a row)."""
    by_name = {member.name: member for member in members}
    generation = {member.name: 0 for member in members}
    children = {name: [] for name in by_name}
    for member in members:
        for parent in (member.father, member.mother):
            if parent in children:
                children[parent].append(member.name)
    for _ in range(len(members) + 1):
        changed = False
        for member in members:
            value = generation[member.name]
            for parent in (member.father, member.mother):
                if parent in generation:
                    value = max(value, generation[parent] + 1)
            if member.father not in generation and member.mother not in generation:
                for child in children[member.name]:
                    value = max(value, generation[child] - 1)
            if value != generation[member.name]:
                generation[member.name] = value
                changed = True
        if not changed:
            break
    return generation


def layout(members: list[Member]) -> dict[str, tuple[float, int]]:
    """Member -> (x, generation): children centered under their parents in their order, a
    founder next to its partner, nodes of a row at least one unit apart."""
    generation = generations(members)
    by_name = {member.name: member for member in members}
    rows: dict[int, list[Member]] = {}
    for member in members:
        rows.setdefault(generation[member.name], []).append(member)
    x: dict[str, float] = {}
    for row in sorted(rows):
        pending = list(rows[row])
        # Children of each couple, under the couple's center
        couples: dict[tuple[str, str], list[Member]] = {}
        for member in pending:
            if member.father in x and member.mother in x:
                couples.setdefault((member.father, member.mother), []).append(member)
        for (father, mother), kids in sorted(
            couples.items(), key=lambda item: (x[item[0][0]] + x[item[0][1]]) / 2
        ):
            center = (x[father] + x[mother]) / 2
            for index, kid in enumerate(kids):
                x[kid.name] = center - (len(kids) - 1) / 2 + index
        # Founders next to their partner, the father to the left of the mother
        for member in pending:
            if member.name in x:
                continue
            partners = [
                (other, member.name == other.father)
                for other in members
                if member.name in (other.father, other.mother)
                for other in [
                    by_name[other.father if member.name == other.mother else other.mother]
                ]
                if other.name in x
            ]
            if partners:
                partner, is_father = partners[0]
                x[member.name] = x[partner.name] + (-1 if is_father else 1)
            else:
                x[member.name] = (
                    (max(x[m.name] for m in pending if m.name in x) + 1)
                    if any(m.name in x for m in pending)
                    else 0.0
                )
        # No two nodes of a row closer than one unit
        ordered = sorted(pending, key=lambda member: x[member.name])
        for previous, member in zip(ordered, ordered[1:]):
            if x[member.name] < x[previous.name] + 1:
                x[member.name] = x[previous.name] + 1
    offset = min(x.values()) if x else 0
    return {name: (value - offset, generation[name]) for name, value in x.items()}


def nodes_and_edges(members: list[Member]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The records of the layout: nodes with their glyph and position, edges between partners
    and from each couple's center to its children."""
    positions = layout(members)
    nodes = [
        {
            "name": member.name,
            "label": member.label
            if len(member.label) <= LABEL_LENGTH
            else member.label[: LABEL_LENGTH - 1] + "…",
            "full_label": member.label,
            "sample": member.sample,
            "x": positions[member.name][0],
            "y": -positions[member.name][1],
            "sex": member.sex if member.sex in ("male", "female") else "unknown",
            "status": {True: "affected", False: "unaffected"}.get(member.affected, "unknown"),
            "observed": member.sample is not None,
        }
        for member in members
    ]
    edges = []
    couples = set()
    for member in members:
        if member.father in positions and member.mother in positions:
            couples.add((member.father, member.mother))
            (xf, yf), (xm, _) = positions[member.father], positions[member.mother]
            edges.append(
                {
                    "x": (xf + xm) / 2,
                    "y": -yf,
                    "x2": positions[member.name][0],
                    "y2": -positions[member.name][1],
                }
            )
    for father, mother in sorted(couples):
        (xf, yf), (xm, ym) = positions[father], positions[mother]
        edges.append({"x": xf, "y": -yf, "x2": xm, "y2": -ym})
    return nodes, edges


DATA_FUNCTION = """function(value, row) {
  const nodes = %(nodes)s;
  const edges = %(edges)s;
  const num = (v) => {
    if (typeof v === "number") return v;
    if (typeof v !== "string" || v === "") return null;
    const m = v.startsWith("<") ? v.match(/data-value="([^"]*)"/) : null;
    const n = parseFloat(m ? m[1] : v);
    return Number.isNaN(n) ? null : n;
  };
  const text = (v) => (typeof v === "string" && v.startsWith("<")) ? (v.match(/data-value="([^"]*)"/) || ["", ""])[1] : (v == null ? "" : String(v));
  const interval = (afd) => {
    const points = [];
    for (const entry of String(afd || "").split(",")) {
      const [af, phred] = entry.split("=");
      if (phred !== undefined) points.push({af: parseFloat(af), w: Math.pow(10, -parseFloat(phred) / 10)});
    }
    const total = points.reduce((s, p) => s + p.w, 0);
    if (!total) return [null, null];
    let cum = 0, lo = null, hi = null;
    for (const p of points) {
      cum += p.w / total;
      if (lo === null && cum >= 0.05) lo = p.af;
      if (hi === null && cum >= 0.95) { hi = p.af; break; }
    }
    return [lo, hi === null ? points[points.length - 1].af : hi];
  };
  const records = edges.map((e) => Object.assign({kind: "edge"}, e));
  for (const node of nodes) {
    const r = Object.assign({kind: "node", text: "no data"}, node);
    if (node.sample) {
      const af = num(row[node.sample + ": AF"]);
      const dp = num(row[node.sample + ": DP"]);
      const gt = text(row[node.sample + ": GT"]);
      const [lo, hi] = interval(text(row[node.sample + ": AFD"]));
      r.af = af;
      r.bar_x = node.x - 0.4;
      r.bar_end = node.x + 0.4;
      r.bar_x2 = af === null ? node.x - 0.4 : node.x - 0.4 + 0.8 * af;
      r.lo_x = lo === null ? null : node.x - 0.4 + 0.8 * lo;
      r.hi_x = hi === null ? null : node.x - 0.4 + 0.8 * hi;
      r.text = (af === null ? "AF ?" : "AF " + af.toFixed(2)) + (dp === null ? "" : " · DP " + dp) + (gt ? " · " + gt : "");
    }
    records.push(r);
  }
  const event = text(row["event"]);
  const p = num(row["P(" + event + ")"]);
  records.push({kind: "title", x: %(center)s, y: 0.75, text: event + (p === null ? "" : " · P = " + p.toFixed(3))});
  return records;
}"""


def plot(members: list[Member]) -> dict[str, Any]:
    """The custom plot of the pedigree: the vega-lite spec and the data function that fills
    the layout with the row's values."""
    nodes, edges = nodes_and_edges(members)
    xs = [node["x"] for node in nodes]
    ys = [node["y"] for node in nodes]
    x_min, x_max = min(xs) - 0.5, max(xs) + 0.5
    y_min, y_max = min(ys) - 0.9, max(ys) + 0.95
    layers = [
        {
            "transform": [{"filter": "datum.kind == 'edge'"}],
            "mark": {"type": "rule", "stroke": MUTED, "strokeWidth": 1},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "x2": {"field": "x2"},
                "y2": {"field": "y2"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "rect", "fill": "#e1e0d9", "height": 5},
            "encoding": {
                "x": {"field": "bar_x", "type": "quantitative"},
                "x2": {"field": "bar_end"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 24},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node' && isValid(datum.af)"}],
            "mark": {"type": "rect", "fill": INK, "height": 5},
            "encoding": {
                "x": {"field": "bar_x", "type": "quantitative"},
                "x2": {"field": "bar_x2"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 24},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node' && isValid(datum.lo_x)"}],
            "mark": {"type": "rule", "stroke": "#0d366b", "strokeWidth": 2},
            "encoding": {
                "x": {"field": "lo_x", "type": "quantitative"},
                "x2": {"field": "hi_x"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 30},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "point", "size": 420, "strokeWidth": 1.5, "tooltip": True},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "shape": {
                    "field": "sex",
                    "type": "nominal",
                    "scale": {
                        "domain": ["male", "female", "unknown"],
                        "range": ["square", "circle", "diamond"],
                    },
                    "legend": None,
                },
                "fill": {
                    "field": "status",
                    "type": "nominal",
                    "scale": {
                        "domain": ["affected", "unaffected", "unknown"],
                        "range": [TEXT, "#ffffff", "#f0efec"],
                    },
                    "legend": None,
                },
                "stroke": {"value": TEXT},
                "strokeDash": {
                    "field": "observed",
                    "type": "nominal",
                    "scale": {"domain": [True, False], "range": [[1, 0], [3, 2]]},
                    "legend": None,
                },
                "tooltip": [
                    {"field": "full_label", "type": "nominal", "title": "member"},
                    {"field": "sample", "type": "nominal"},
                    {"field": "text", "type": "nominal", "title": "call"},
                ],
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "text", "fontSize": 10, "color": TEXT, "yOffset": 44},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "text": {"field": "label"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "text", "fontSize": 9, "color": TEXT, "yOffset": 56},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "text": {"field": "text"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'title'"}],
            "mark": {"type": "text", "fontSize": 11, "color": TEXT, "fontWeight": "bold"},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "text": {"field": "text"},
            },
        },
    ]
    spec = {
        "$schema": VEGA_LITE,
        "width": int((x_max - x_min) * X_STEP),
        "height": int((y_max - y_min) * Y_STEP),
        "encoding": {
            "x": {"scale": {"domain": [x_min, x_max]}, "axis": None},
            "y": {"scale": {"domain": [y_min, y_max]}, "axis": None},
        },
        "layer": layers,
        "config": CELL_CONFIG,
    }
    data = DATA_FUNCTION % {
        "nodes": json.dumps(nodes),
        "edges": json.dumps(edges),
        "center": (x_min + x_max) / 2,
    }
    return {"data": data, "spec": spec}


def column(members: list[Member]) -> Column:
    """The ``pedigree`` column: the plot in the row details."""
    return Column(custom_plot=plot(members), display="detail")
