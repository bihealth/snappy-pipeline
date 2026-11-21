# -*- coding: utf-8 -*-
"""The pedigree of a call as a plot in the row details: who carries the variant, and how surely.

The members of the entity's pedigree, from the sample sheet, are laid out once per table:
generations as rows, children below their parents, partners side by side, with the unsequenced donors that connect members and the unknown parents of members with one parent drawn hollow.
Per call, the plot reads the per-sample columns of the row: the allele frequency as a bar under the node with the 90 % credible interval from the allele frequency density, below it the reads as a diverging bar (for the alternative allele to the right in blue, for the reference to the left in gray, darker for stronger evidence, pale for uncertain mapping, scaled to the member with the most reads, on a gray background of that width), the depth and the genotype as text, and the call's event with its posterior as the title.
The glyphs follow the pedigree convention: square for male, circle for female, diamond for unknown sex, filled when affected;
    the stroke carries the sex color (aqua, orange) and the fill of an affected member repeats it; a legend sits under the pedigree.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

from snappy_wrappers.tools.qc_summary import (
    CELL_CONFIG,
    COLORS,
    INK,
    MUTED,
    TEXT,
    VEGA_LITE,
    Column,
)

#: Pixels per unit of the layout, sideways and downwards
X_STEP, Y_STEP = 96, 124

#: Characters of a member's label under its node; longer labels are cut, the tooltip has them
LABEL_LENGTH = 14

#: The stroke of a node by sex: aqua for male, orange for female (the shape says it as well),
#: neither the blue of the evidence nor the usual blue and pink
SEX_COLORS = {"male": COLORS[2], "female": COLORS[1], "unknown": MUTED}

#: The legend under the pedigree: (shape, sex, status, observed, label), in two rows
LEGEND = (
    ("male", "unknown", True, "male"),
    ("female", "unknown", True, "female"),
    ("unknown", "unknown", True, "unknown sex"),
    ("male", "affected", True, "affected"),
    ("male", "unaffected", True, "unaffected"),
    ("male", "unknown", False, "no data"),
)

#: Width in rows that the legend needs
LEGEND_WIDTH = 3.3

#: The fill of a node by affected status: the sex color when affected, white when unaffected,
#: a light gray when unknown
UNAFFECTED_FILL, UNKNOWN_FILL = "#ffffff", "#f0efec"

#: The reads of a member by strength of evidence (equal, barely, positive, strong, very
#: strong): for the alternative allele on the blue ramp, for the reference on a gray ramp
ALT_RAMP = ("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#0d366b")
REF_RAMP = ("#f0efec", "#e1e0d9", "#c3c2b7", "#898781", "#52514e")

#: How far below a couple the sibship line runs (in rows): just above the children, so that
#: the drop from the couple's center is the one long line and the children hang on stubs
SIBSHIP = 0.82

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


def couples_of(members: list[Member]) -> dict[str, list[str]]:
    """Member -> its partners (the other parent of each of its children), in input order."""
    partners: dict[str, list[str]] = {member.name: [] for member in members}
    for member in members:
        father, mother = member.father, member.mother
        if father in partners and mother in partners:
            if mother not in partners[father]:
                partners[father].append(mother)
            if father not in partners[mother]:
                partners[mother].append(father)
    return partners


def layout(members: list[Member]) -> dict[str, tuple[float, int]]:
    """Member -> (x, generation): children centered under their parents in their order, a
    couple side by side with the father to the left, and no node on another couple's line.

    A row is laid out in blocks: a couple of the row is one block, any other member its own;
    a child's block starts under the center of its parents, a founder's next to its partner,
    and blocks that overlap move to the right in the order of their positions."""
    generation = generations(members)
    by_name = {member.name: member for member in members}
    partners = couples_of(members)
    rows: dict[int, list[Member]] = {}
    for member in members:
        rows.setdefault(generation[member.name], []).append(member)
    x: dict[str, float] = {}
    for row in sorted(rows):
        names = [member.name for member in rows[row]]
        in_row = set(names)
        # Blocks: a couple of this row (father, mother), else one member
        blocks: list[list[str]] = []
        taken: set[str] = set()
        for name in names:
            if name in taken:
                continue
            partner = next((p for p in partners[name] if p in in_row and p not in taken), None)
            if partner is None:
                blocks.append([name])
                taken.add(name)
                continue
            father = (
                name
                if any(by_name[c].father == name for c in by_name if by_name[c].mother == partner)
                else partner
            )
            mother = partner if father == name else name
            blocks.append([father, mother])
            taken.update((father, mother))
        # The start of a block: under the parents of its members, else after the last block
        starts: list[float | None] = []
        for block in blocks:
            centers = []
            for name in block:
                member = by_name[name]
                if member.father in x and member.mother in x:
                    centers.append((x[member.father] + x[member.mother]) / 2)
            if centers:
                # Children of one couple spread around its center, in their order
                siblings = (
                    [
                        c
                        for c in names
                        if by_name[c].father == by_name[block[0]].father
                        and by_name[c].mother == by_name[block[0]].mother
                    ]
                    if by_name[block[0]].father in x
                    else [block[0]]
                )
                index = siblings.index(block[0]) if block[0] in siblings else 0
                starts.append(centers[0] - (len(siblings) - 1) / 2 + index)
            else:
                starts.append(None)
        # Children first (their start is known), then the others after the rightmost block
        order = sorted(range(len(blocks)), key=lambda i: (starts[i] is None, starts[i] or 0))
        right = None
        for i in order:
            if starts[i] is None:
                starts[i] = 0.0 if right is None else right + 1
            right = (
                starts[i] + len(blocks[i]) - 1
                if right is None
                else max(right, starts[i] + len(blocks[i]) - 1)
            )
        # Blocks never overlap: shift to the right in the order of their starts
        placed = sorted(range(len(blocks)), key=lambda i: starts[i])
        edge = None
        for i in placed:
            if edge is not None and starts[i] < edge + 1:
                starts[i] = edge + 1
            edge = starts[i] + len(blocks[i]) - 1
            for offset, name in enumerate(blocks[i]):
                x[name] = starts[i] + offset
    offset = min(x.values()) if x else 0
    return {name: (value - offset, generation[name]) for name, value in x.items()}


def nodes_and_edges(members: list[Member]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The records of the layout: nodes with their glyph and position, and the line segments
    between partners, down to the sibship line and down to each child."""
    positions = layout(members)
    nodes = []
    for member in members:
        sex = member.sex if member.sex in ("male", "female") else "unknown"
        status = {True: "affected", False: "unaffected"}.get(member.affected, "unknown")
        fills = {
            "affected": SEX_COLORS[sex],
            "unaffected": UNAFFECTED_FILL,
            "unknown": UNKNOWN_FILL,
        }
        nodes.append(
            {
                "name": member.name,
                "label": member.label
                if len(member.label) <= LABEL_LENGTH
                else member.label[: LABEL_LENGTH - 1] + "…",
                "full_label": member.label,
                "sample": member.sample,
                "x": positions[member.name][0],
                "y": -positions[member.name][1],
                "sex": sex,
                "status": status,
                "stroke_color": SEX_COLORS[sex],
                "fill_color": fills[status],
                "observed": member.sample is not None,
            }
        )
    # Orthogonal lines, as pedigrees are drawn: the partners joined, a drop from the middle of
    # their line to the sibship line half-way down, and a drop from it to each child
    couples: dict[tuple[str, str], list[str]] = {}
    for member in members:
        if member.father in positions and member.mother in positions:
            couples.setdefault((member.father, member.mother), []).append(member.name)
    edges = []
    for (father, mother), children in sorted(couples.items()):
        (xf, yf), (xm, _) = positions[father], positions[mother]
        top, center = -yf, (xf + xm) / 2
        sibship = top - SIBSHIP
        edges.append({"x": xf, "y": top, "x2": xm, "y2": top})
        edges.append({"x": center, "y": top, "x2": center, "y2": sibship})
        xs = [positions[child][0] for child in children] + [center]
        if max(xs) > min(xs):
            edges.append({"x": min(xs), "y": sibship, "x2": max(xs), "y2": sibship})
        for child in children:
            x, y = positions[child]
            edges.append({"x": x, "y": sibship, "x2": x, "y2": -y})
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
  const strengths = {E: 0, B: 1, P: 2, S: 3, V: 4};
  const strengthNames = ["equal", "barely", "positive", "strong", "very strong"];
  const altRamp = %(alt_ramp)s;
  const refRamp = %(ref_ramp)s;
  const reads = (textValue) => {
    const sides = String(textValue || "").split("|");
    const parse = (side) => {
      const out = [];
      for (const m of (side || "").matchAll(/(\\d+)([A-Za-z])/g)) {
        const code = m[2].toUpperCase();
        if (code in strengths) out.push({strength: strengths[code], count: parseInt(m[1]), certain: m[2] === code});
      }
      return out.sort((a, b) => b.strength - a.strength);
    };
    return {alt: parse(sides[0]), ref: parse(sides[1])};
  };
  const observations = {};
  let most = 0;
  for (const node of nodes) {
    if (!node.sample) continue;
    observations[node.name] = reads(text(row[node.sample + ": observations"]));
    for (const side of ["alt", "ref"]) {
      most = Math.max(most, observations[node.name][side].reduce((s, o) => s + o.count, 0));
    }
  }
  const unit = most > 0 ? 0.4 / most : 0;
  for (const node of nodes) {
    const r = Object.assign({kind: "node", text: "no data"}, node);
    if (node.sample) {
      const obs = observations[node.name];
      for (const [side, sign, ramp] of [["alt", 1, altRamp], ["ref", -1, refRamp]]) {
        let at = node.x;
        for (const o of obs[side]) {
          const end = at + sign * o.count * unit;
          records.push({kind: "obs", x: Math.min(at, end), x2: Math.max(at, end), y: node.y,
            fill: ramp[o.strength], opacity: o.certain ? 1 : 0.45,
            label: `${side === "alt" ? "alternative allele" : "reference"} · ${strengthNames[o.strength]} · ${o.count} reads${o.certain ? "" : " (uncertain mapping)"}`});
          at = end;
        }
      }
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


def legend_records(x_min: float, y_min: float) -> list[dict[str, Any]]:
    """The glyphs of the legend with their labels: sex in the first row, status in the second,
    and the reading of the bars as a line of text."""
    records = []
    for index, (sex, status, observed, label) in enumerate(LEGEND):
        row, column = divmod(index, 3)
        fills = {
            "affected": SEX_COLORS[sex],
            "unaffected": UNAFFECTED_FILL,
            "unknown": UNKNOWN_FILL,
        }
        records.append(
            {
                "kind": "legend",
                "x": x_min + 0.25 + column * 1.05,
                "y": y_min - 0.25 - row * 0.3,
                "sex": sex,
                "stroke_color": SEX_COLORS[sex]
                if status == "unknown" or sex != "male"
                else SEX_COLORS["male"],
                "fill_color": fills[status],
                "observed": observed,
                "text": label,
            }
        )
    records.append(
        {
            "kind": "note",
            "x": x_min + 0.1,
            "y": y_min - 0.9,
            "text": "bars: allele frequency with its 90 % interval; reads ← reference · alternative →, darker = stronger evidence",
        }
    )
    return records


def plot(members: list[Member]) -> dict[str, Any]:
    """The custom plot of the pedigree: the vega-lite spec and the data function that fills
    the layout with the row's values."""
    nodes, edges = nodes_and_edges(members)
    xs = [node["x"] for node in nodes]
    ys = [node["y"] for node in nodes]
    x_min, x_max = min(xs) - 0.5, max(max(xs) + 0.5, min(xs) - 0.5 + LEGEND_WIDTH)
    y_min, y_max = min(ys) - 0.9, max(ys) + 0.95
    legend = legend_records(x_min, y_min)
    y_min -= 1.1
    layers = [
        {
            "transform": [{"filter": "datum.kind == 'edge'"}],
            "mark": {"type": "rule", "stroke": MUTED, "strokeWidth": 1, "opacity": 1},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "x2": {"field": "x2"},
                "y2": {"field": "y2"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "rect", "fill": "#e1e0d9", "height": 5, "opacity": 1},
            "encoding": {
                "x": {"field": "bar_x", "type": "quantitative"},
                "x2": {"field": "bar_end"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 20},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node' && isValid(datum.af)"}],
            "mark": {"type": "rect", "fill": INK, "height": 5, "opacity": 1},
            "encoding": {
                "x": {"field": "bar_x", "type": "quantitative"},
                "x2": {"field": "bar_x2"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 20},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node' && isValid(datum.lo_x)"}],
            "mark": {"type": "rule", "stroke": "#0d366b", "strokeWidth": 2, "opacity": 1},
            "encoding": {
                "x": {"field": "lo_x", "type": "quantitative"},
                "x2": {"field": "hi_x"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 26},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node' && datum.sample"}],
            "mark": {"type": "rect", "fill": "#e1e0d9", "height": 6, "opacity": 1},
            "encoding": {
                "x": {"field": "bar_x", "type": "quantitative"},
                "x2": {"field": "bar_end"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 35},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'obs'"}],
            "mark": {"type": "rect", "height": 6, "tooltip": True, "opacity": 1},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "x2": {"field": "x2"},
                "y": {"field": "y", "type": "quantitative"},
                "yOffset": {"value": 35},
                "fill": {"field": "fill", "type": "nominal", "scale": None},
                "fillOpacity": {"field": "opacity", "type": "quantitative", "scale": None},
                "tooltip": [{"field": "label", "type": "nominal", "title": "reads"}],
            },
        },
        {
            # The center of the reads; a rule needs its y2, or it runs to the zero line
            "transform": [{"filter": "datum.kind == 'node' && datum.sample"}],
            "mark": {
                "type": "rule",
                "stroke": TEXT,
                "strokeWidth": 1,
                "opacity": 1,
                "yOffset": 31,
                "y2Offset": 39,
            },
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "y2": {"field": "y"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {
                "type": "point",
                "size": 420,
                "strokeWidth": 1.5,
                "tooltip": True,
                "opacity": 1,
            },
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
                "fill": {"field": "fill_color", "type": "nominal", "scale": None},
                "stroke": {"field": "stroke_color", "type": "nominal", "scale": None},
                "strokeDash": {
                    "field": "observed",
                    "type": "nominal",
                    "scale": {"domain": [True, False], "range": [[1, 0], [3, 2]]},
                    "legend": None,
                },
                "tooltip": [
                    {"field": "full_label", "type": "nominal", "title": "member"},
                    {"field": "sex", "type": "nominal"},
                    {"field": "status", "type": "nominal"},
                    {"field": "sample", "type": "nominal"},
                    {"field": "text", "type": "nominal", "title": "call"},
                ],
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "text", "fontSize": 10, "color": TEXT, "yOffset": 50},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "text": {"field": "label"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'node'"}],
            "mark": {"type": "text", "fontSize": 9, "color": TEXT, "yOffset": 62},
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
        {
            "transform": [{"filter": "datum.kind == 'legend'"}],
            "mark": {"type": "point", "size": 110, "strokeWidth": 1.5, "opacity": 1},
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
                "fill": {"field": "fill_color", "type": "nominal", "scale": None},
                "stroke": {"field": "stroke_color", "type": "nominal", "scale": None},
                "strokeDash": {
                    "field": "observed",
                    "type": "nominal",
                    "scale": {"domain": [True, False], "range": [[1, 0], [3, 2]]},
                    "legend": None,
                },
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'legend'"}],
            "mark": {"type": "text", "fontSize": 9, "color": TEXT, "align": "left", "dx": 10},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "text": {"field": "text"},
            },
        },
        {
            "transform": [{"filter": "datum.kind == 'note'"}],
            "mark": {"type": "text", "fontSize": 8, "color": MUTED, "align": "left"},
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
        "edges": json.dumps(edges + legend),
        "center": (x_min + x_max) / 2,
        "alt_ramp": json.dumps(list(ALT_RAMP)),
        "ref_ramp": json.dumps(list(REF_RAMP)),
    }
    return {"data": data, "spec": spec}


def column(members: list[Member]) -> Column:
    """The ``pedigree`` column: the plot in the row details."""
    return Column(custom_plot=plot(members), display="detail")
