# -*- coding: utf-8 -*-
"""The views of the calls site.

With varlociraptor calls,
one view per event of the scenario (not artifact and absent)
holds the calls whose most probable event it is, ranked by that posterior,
then by gnomAD frequency (unobserved first), then by impact;
the site opens with the first event view that has calls.
The per-member events ``denovo_<member>`` share the color of their kind.
Annotated calls without events get a ``ranked`` view by impact, frequency and quality.
``all`` lists every call in genomic order, without the observation and density plots,
the other annotations and the inspect link, which the event views and the result TSV keep;
a row links to the same call in the other views.
"""

from __future__ import annotations

import re
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import COLORS, MUTED, Column, Summary, Table, pills

from . import annotation, varlociraptor
from .vcf import Vcf

#: The events of the preset scenarios in the fixed categorical order; artifact and absent muted
EVENT_COLORS = {
    "present": COLORS[0],
    "denovo": COLORS[1],
    "recessive": COLORS[2],
    "dominant": COLORS[3],
    "somatic_tumor_high": COLORS[1],
    "somatic_tumor_medium": COLORS[3],
    "somatic_tumor_low": COLORS[4],
    "artifact": MUTED,
    "absent": MUTED,
}

#: Columns that only the event views and the result TSV carry:
#: the per-sample plots, the other annotations and the inspect link
HEAVY = ("other annotations", "inspect", "inspect_url")

#: The link from a row to the same call in another view, anchored on the position
LINK_COLUMN = "position"


def event_names(events: list[str]) -> list[str]:
    return [event.lower() for event in events]


def event_kind(event: str) -> str:
    """The kind of an event, for its color:
    ``denovo`` for the germline preset's per-member events ``denovo_<member>``,
    else the event itself."""
    return "denovo" if event.startswith("denovo_") else event


def event_label(event: str) -> str:
    """The event with the member it names, for the germline preset's ``denovo_<member>``:
    ``denovo_idx (de novo in idx)``."""
    if event.startswith("denovo_"):
        return f"{event} (de novo in {event.removeprefix('denovo_')})"
    return event


def event_colors(events: list[str]) -> tuple[list[str], list[str]]:
    """The pills domain and colors of the events:
    the preset events and kinds keep their color, others take the next free categorical slot."""
    domain = event_names(events)
    used = {
        EVENT_COLORS[event_kind(event)] for event in domain if event_kind(event) in EVENT_COLORS
    }
    free = [color for color in COLORS if color not in used]
    colors = []
    for event in domain:
        if event_kind(event) in EVENT_COLORS:
            colors.append(EVENT_COLORS[event_kind(event)])
        else:
            colors.append(free.pop(0) if free else MUTED)
    return domain, colors


def event_of(row: dict[str, Any], events: list[str]) -> str:
    """The most probable event of a row (ties: the first in header order); empty without any."""
    best, best_probability = "", -1.0
    for event in events:
        probability = row.get(varlociraptor.probability_column(event))
        if probability is not None and probability > best_probability:
            best, best_probability = event.lower(), probability
    return best


def event_column(events: list[str]) -> Column:
    domain, colors = event_colors(events)
    return pills(domain, colors)


def _natural(name: str) -> tuple:
    return tuple(int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name))


def contig_ranks(frame: pl.DataFrame, contigs: list[str]) -> dict[str, int]:
    """Chromosome -> rank: the header's contig order, contigs it lacks after them in natural order."""
    known = {contig: index for index, contig in enumerate(contigs)}
    others = sorted(
        (c for c in frame["chromosome"].unique().to_list() if c not in known and c is not None),
        key=_natural,
    )
    return known | {contig: len(known) + index for index, contig in enumerate(others)}


def _with_ranks(frame: pl.DataFrame, contigs: list[str]) -> pl.DataFrame:
    ranks = contig_ranks(frame, contigs)
    frame = frame.with_columns(
        pl.col("chromosome")
        .replace_strict(ranks, default=len(ranks), return_dtype=pl.Int64)
        .alias("_contig")
    )
    if "impact" in frame.columns:
        impacts = {impact: rank for rank, impact in enumerate(annotation.IMPACTS)}
        frame = frame.with_columns(
            pl.col("impact")
            .replace_strict(impacts, default=len(impacts), return_dtype=pl.Int64)
            .alias("_impact")
        )
    return frame


def genomic_order(frame: pl.DataFrame, contigs: list[str]) -> pl.DataFrame:
    if frame.height == 0:
        return frame
    return (
        _with_ranks(frame, contigs)
        .sort("_contig", "position")
        .drop("_contig", "_impact", strict=False)
    )


def ranked(frame: pl.DataFrame, contigs: list[str], first: str | None = None) -> pl.DataFrame:
    """``frame`` by ``first`` (a probability, descending), then gnomAD frequency (unobserved first),
    impact, quality (descending) and position."""
    if frame.height == 0:
        return frame
    frame = _with_ranks(frame, contigs)
    keys: list[tuple[str, bool, bool]] = []  # column, descending, nulls last
    if first:
        keys.append((first, True, True))
    if "gnomAD AF" in frame.columns:
        keys.append(("gnomAD AF", False, False))
    if "_impact" in frame.columns:
        keys.append(("_impact", False, False))
    if "quality" in frame.columns:
        keys.append(("quality", True, True))
    keys += [("_contig", False, False), ("position", False, False)]
    return frame.sort(
        [key for key, _, _ in keys],
        descending=[descending for _, descending, _ in keys],
        nulls_last=[nulls_last for _, _, nulls_last in keys],
    ).drop("_contig", "_impact", strict=False)


def _ranking_text(frame: pl.DataFrame, posterior: bool) -> str:
    parts = ["that posterior"] if posterior else []
    if "gnomAD AF" in frame.columns:
        parts.append("gnomAD frequency (unobserved first)")
    if "impact" in frame.columns:
        parts.append("impact")
    if not posterior and "quality" in frame.columns:
        parts.append("quality")
    parts.append("position")
    return ", then ".join(parts)


def build(
    vcf: Vcf,
    frame: pl.DataFrame,
    columns: dict[str, Column],
    page_size: int,
    max_in_memory_rows: int | None,
    description,
) -> Summary:
    """The summary of the site: ``all`` as its table, the event or ranked views, the default."""
    events = varlociraptor.events(vcf)
    for_events = [e for e in events if e.upper() not in varlociraptor.AGAINST_EVENTS]
    heavy = tuple(
        name
        for name, column in columns.items()
        if name in frame.columns and (column.custom_plot or name in HEAVY)
    )
    settings = {
        "columns": columns,
        "pinned": ("chromosome", "position"),
        "page_size": page_size,
        "max_in_memory_rows": max_in_memory_rows,
    }
    views: list[Table] = []
    for event in event_names(for_events):
        calls = frame.filter(pl.col("event") == event) if frame.height else frame
        calls = ranked(calls, vcf.contigs, varlociraptor.probability_column(event))
        lead = (
            f"The {calls.height} calls whose most probable event is {event_label(event)}, "
            f"ranked by {_ranking_text(calls, True)}."
        )
        views.append(
            Table(
                event,
                f"{event} calls",
                calls,
                desc=description.with_lead(
                    lead, navigation=["A row links to the same call in all calls."]
                ).markdown(),
                links={"all calls": (LINK_COLUMN, f"all/{LINK_COLUMN}", False)},
                **settings,
            )
        )
    if not events and any(name in frame.columns for name in ("impact", "gnomAD AF", "quality")):
        calls = ranked(frame, vcf.contigs)
        views.append(
            Table(
                "ranked",
                "ranked calls",
                calls,
                desc=description.with_lead(
                    f"All calls ranked by {_ranking_text(calls, False)}.",
                    navigation=["A row links to the same call in all calls."],
                ).markdown(),
                links={"all calls": (LINK_COLUMN, f"all/{LINK_COLUMN}", False)},
                **settings,
            )
        )
    kinds = []
    if any(columns[name].custom_plot for name in heavy):
        kinds.append("the observation and density plots")
    if "other annotations" in heavy:
        kinds.append("the other annotations")
    if "inspect" in heavy:
        kinds.append("the inspect link")
    navigation = []
    if kinds:
        listed = ", ".join(kinds[:-1]) + (" and " if len(kinds) > 1 else "") + kinds[-1]
        navigation.append(f"{listed.capitalize()} are in the other views and in the result TSV.")
    if views:
        navigation.append("A row links to the same call in its other views.")
    if max_in_memory_rows is not None and frame.height > max_in_memory_rows:
        navigation.append(
            f"With {frame.height} rows, above the in-memory limit of {max_in_memory_rows}, "
            "this view is paged without sorting; the other views sort."
        )
    table = Table(
        "all",
        "all calls",
        genomic_order(frame, vcf.contigs),
        desc=description.with_lead("All calls in genomic order.", navigation=navigation).markdown(),
        links={view.name: (LINK_COLUMN, f"{view.name}/{LINK_COLUMN}", True) for view in views},
        exclude=heavy,
        **settings,
    )
    default = next((view.name for view in views if view.frame.height), None)
    return Summary(table=table, views=views, default_view=default)
