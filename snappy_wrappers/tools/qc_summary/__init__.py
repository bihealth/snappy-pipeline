# -*- coding: utf-8 -*-
"""Tables and charts of the QC files of a task's libraries, rendered with datavzrd.

A format module (``samtools``, ``alfred``, ``picard``, ``bcftools_stats``) parses the files of one
tool for all libraries into a :class:`Summary`: the summary table, one row per library, which the
task writes as its result TSV, and the views of the site, tables and charts.
:func:`render` writes the TSV and runs datavzrd, which turns the views into static HTML that the
Snakemake report embeds (``snappy report``).

Frames are polars frames with one row per library and observation; the ``library`` column comes
first, so the site pins it. The QC files are small (kilobytes per library), so the frames are
built eagerly. Distributions are thinned to :data:`MAX_POINTS` per library before plotting.
"""

from __future__ import annotations

import importlib
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass, field

import polars as pl

#: Points per library that a chart keeps of a distribution
MAX_POINTS = 120

#: Series that the legend of a chart lists; the tooltip names the rest
LEGEND_SYMBOLS = 20


@dataclass(frozen=True)
class Column:
    """Rendering of one table column: its label, precision, cell plot and visibility."""

    label: str | None = None
    precision: int | None = None
    #: A datavzrd cell plot, for example ``{"heatmap": {...}}`` or ``{"bars": {...}}``
    plot: dict | None = None
    #: Shown in the details of a row only
    detail: bool = False

    def spec(self) -> dict:
        spec: dict = {}
        if self.label:
            spec["label"] = self.label
        if self.precision is not None:
            spec["precision"] = self.precision
        if self.plot:
            spec["plot"] = self.plot
        if self.detail:
            spec["display-mode"] = "detail"
        return spec


def fraction(label: str | None = None, detail: bool = False, high: float = 1.0) -> Column:
    """A fraction or percentage from 0 to ``high``, shaded by its value."""
    # datavzrd takes the domain of a heatmap as text (its schema), also for a linear scale
    domain = ["0", str(high)]
    plot = {"heatmap": {"scale": "linear", "domain": domain, "color-scheme": "blues"}}
    return Column(label=label, precision=3 if high == 1.0 else 2, plot=plot, detail=detail)


def bars(label: str | None = None, detail: bool = False, precision: int = 0) -> Column:
    """A count or measure, with a bar proportional to its value."""
    return Column(
        label=label, precision=precision, plot={"bars": {"scale": "linear"}}, detail=detail
    )


@dataclass
class Table:
    """A table of the site; ``columns`` sets the rendering of the columns it names."""

    name: str
    title: str
    frame: pl.DataFrame
    desc: str = ""
    columns: dict[str, Column] = field(default_factory=dict)
    #: Columns kept at the left when the table scrolls
    pinned: tuple[str, ...] = ("library",)
    #: Rows per page
    page_size: int = 25


@dataclass
class Chart:
    """A chart of the site from a long frame: one line or bar series per ``color`` value."""

    name: str
    title: str
    frame: pl.DataFrame
    x: str
    y: str
    desc: str = ""
    color: str = "library"
    mark: str = "line"
    #: ``quantitative`` or ``nominal`` (categories on the x axis)
    x_type: str = "quantitative"
    #: Column whose values become side-by-side panels, for example the read of a pair
    facet: str | None = None
    log_y: bool = False
    x_title: str | None = None
    y_title: str | None = None

    def spec(self) -> dict:
        """The vega-lite specification of the chart."""
        x: dict = {"field": self.x, "type": self.x_type, "title": self.x_title or self.x}
        y: dict = {"field": self.y, "type": "quantitative", "title": self.y_title or self.y}
        if self.log_y:
            y["scale"] = {"type": "log"}
        encoding: dict = {
            "x": x,
            "y": y,
            "color": {
                "field": self.color,
                "type": "nominal",
                "legend": {"symbolLimit": LEGEND_SYMBOLS},
            },
            "tooltip": [
                {"field": self.color, "type": "nominal"},
                {"field": self.x, "type": self.x_type},
                {"field": self.y, "type": "quantitative"},
            ],
        }
        if self.mark == "bar" and self.x_type == "nominal":
            encoding["xOffset"] = {"field": self.color}
        if self.facet:
            encoding["column"] = {"field": self.facet, "type": "nominal", "title": None}
        return {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "title": self.title,
            "mark": {"type": self.mark, "tooltip": True},
            "encoding": encoding,
            "width": 320 if self.facet else 640,
            "height": 320,
        }


@dataclass
class Summary:
    """What one format makes of the files: the summary table and the views of the site."""

    table: Table
    views: list[Table | Chart]


def summarize(format: str, files: dict[str, list[str]], labels: list[str]) -> Summary:
    """Parse ``files`` (paths by kind, one per label in the same order) with the format module."""
    module = importlib.import_module(f"{__name__}.{format}")
    return module.summarize(files, labels)


def render(summary: Summary, title: str, tsv: str, site: str) -> None:
    """Write the summary table to ``tsv`` and the datavzrd site of all views to ``site``."""
    os.makedirs(os.path.dirname(os.path.abspath(tsv)), exist_ok=True)
    summary.table.frame.write_csv(tsv, separator="\t")
    with tempfile.TemporaryDirectory() as tmpdir:
        views = [summary.table, *summary.views]
        for view in views:
            view.frame.write_csv(os.path.join(tmpdir, f"{view.name}.tsv"), separator="\t")
        config_path = os.path.join(tmpdir, "config.json")
        with open(config_path, "wt") as f:
            json.dump(config(views, title, tmpdir), f, indent=1)
        subprocess.run(
            ["datavzrd", config_path, "--output", site, "--overwrite-output"], check=True
        )


def config(views: list[Table | Chart], title: str, data_dir: str) -> dict:
    """The datavzrd configuration of ``views``, whose frames are TSV files in ``data_dir``."""
    datasets = {
        view.name: {"path": os.path.join(data_dir, f"{view.name}.tsv"), "separator": "\t"}
        for view in views
    }
    return {
        "name": title,
        "default-view": views[0].name,
        "datasets": datasets,
        "views": {view.name: _view(view) for view in views},
    }


def _view(view: Table | Chart) -> dict:
    spec: dict = {"dataset": view.name, "desc": view.desc or view.title}
    if isinstance(view, Chart):
        spec["render-plot"] = {"spec": json.dumps(view.spec())}
        return spec
    columns = {
        name: {"display-mode": "pinned"} for name in view.pinned if name in view.frame.columns
    }
    for name, column in view.columns.items():
        if name in view.frame.columns:
            columns[name] = columns.get(name, {}) | column.spec()
    spec["page-size"] = view.page_size
    spec["render-table"] = {"columns": columns}
    return spec


def thin(frame: pl.DataFrame, x: str, y: str, by: tuple[str, ...] = ("library",)) -> pl.DataFrame:
    """Bin the distribution ``y`` over ``x`` of each group to at most :data:`MAX_POINTS` points,
    summing ``y`` per bin, so that a chart of many libraries stays small."""
    groups = [name for name in by if name in frame.columns]
    if frame.is_empty():
        return frame
    points = frame.group_by(groups).agg(pl.col(x).n_unique().alias("n"))["n"].max()
    if int(points) <= MAX_POINTS:  # type: ignore[arg-type]
        return frame
    low, high = float(frame[x].min()), float(frame[x].max())  # type: ignore[arg-type]
    width = (high - low) / MAX_POINTS
    return (
        frame.with_columns((((pl.col(x) - low) // width) * width + low).alias(x))
        .group_by([*groups, x], maintain_order=True)
        .agg(pl.col(y).sum())
    )


def concat(frames: list[pl.DataFrame], columns: list[str]) -> pl.DataFrame:
    """Stack the non-empty ``frames`` and order their ``columns``; empty without any."""
    frames = [frame for frame in frames if not frame.is_empty()]
    if not frames:
        return pl.DataFrame(schema={column: pl.Null for column in columns})
    return pl.concat(frames, how="diagonal").select(columns)


def number(text: str) -> float | int:
    """An int or a float from its text."""
    value = float(text)
    return int(value) if value.is_integer() else value


def numeric(frame: pl.DataFrame, columns: list[str] | None = None) -> pl.DataFrame:
    """Cast the text ``columns`` (default: all) that hold only numbers to integers or floats;
    empty cells become null, columns with other text stay as they are."""
    result = frame
    for name in columns or frame.columns:
        column = frame[name]
        if column.dtype != pl.String:
            continue
        blanks = (column == "").sum() + column.null_count()
        for dtype in (pl.Int64, pl.Float64):
            cast = column.cast(dtype, strict=False)
            if cast.null_count() == blanks:
                result = result.with_columns(cast)
                break
    return result


def ratio(numerator, denominator) -> float | None:
    """``numerator / denominator``, ``None`` without both."""
    if numerator is None or denominator is None or not denominator:
        return None
    return float(numerator) / float(denominator)


def from_rows(rows: list[dict], columns: list[str]) -> pl.DataFrame:
    """A frame of ``rows`` (dicts) with ``columns`` first, in order, and the other keys after."""
    if not rows:
        return pl.DataFrame(schema={column: pl.Null for column in columns})
    frame = pl.from_dicts(rows, infer_schema_length=None)
    rest = [column for column in frame.columns if column not in columns]
    return frame.select(
        *(pl.col(c) if c in frame.columns else pl.lit(None).alias(c) for c in columns), *rest
    )
