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

The charts follow Tufte: the data carries the ink. No gridlines, hairline axes, thin lines, no
chart frame; a table cell shows its value as a tick on the range of its column (:func:`ticks`)
or as a bar (:func:`bars`); a distribution is one thin line per library; the libraries of a
summary appear as a strip of ticks per metric (:func:`metric_strips`), so an outlier stands out
and the tooltip names it. Up to :data:`MAX_COLORED` libraries get a color each, in a fixed
order, with a legend; more libraries share one muted ink and are told apart by the tooltip.
"""

from __future__ import annotations

import importlib
import json
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass, field

import polars as pl

#: Points per library that a chart keeps of a distribution
MAX_POINTS = 120

#: Libraries that get a color of their own (and a legend); more share one ink
MAX_COLORED = 6

#: Categorical colors in their fixed order (blue, orange, aqua, yellow, magenta, green), a
#: palette whose adjacent pairs stay apart under color vision deficiencies
COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")

#: The one ink of many libraries, and the inks of text and axes
INK = "#2a78d6"

#: Name of the datavzrd link from a mark of a chart to the row of its series in a table
LINK = "row"
TEXT = "#52514e"
AXIS = "#c3c2b7"

#: Vega-lite configuration of every chart: no grid, hairline axes, no frame, muted text
CHART_CONFIG = {
    "view": {"stroke": None},
    "axis": {
        "grid": False,
        "domainColor": AXIS,
        "tickColor": AXIS,
        "labelColor": TEXT,
        "titleColor": TEXT,
        "titleFontWeight": "normal",
        "labelFontSize": 11,
        "titleFontSize": 11,
    },
    "legend": {
        "labelColor": TEXT,
        "titleColor": TEXT,
        "symbolLimit": MAX_COLORED,
        "labelLimit": 400,
        "title": None,
    },
    "title": {"color": TEXT, "fontWeight": "normal", "fontSize": 13, "anchor": "start"},
    "header": {"labelColor": TEXT, "titleColor": TEXT, "labelFontSize": 11, "titleFontSize": 11},
}


@dataclass(frozen=True)
class Column:
    """Rendering of one table column: its label, precision, cell plot and visibility."""

    label: str | None = None
    precision: int | None = None
    #: A datavzrd cell plot, for example ``{"ticks": {...}}`` or ``{"bars": {...}}``
    plot: dict | None = None
    #: A vega-lite plot per cell: ``{"data": <JS function of value and row>, "spec": {...}}``
    custom_plot: dict | None = None
    #: Shown in the details of a row only
    detail: bool = False
    #: Not shown at all (kept for links, tooltips and the TSV)
    hidden: bool = False

    def spec(self) -> dict:
        spec: dict = {}
        if self.label:
            spec["label"] = self.label
        if self.precision is not None:
            spec["precision"] = self.precision
        if self.plot:
            spec["plot"] = self.plot
        if self.custom_plot:  # datavzrd takes the vega-lite spec of a cell plot as text
            plot = dict(self.custom_plot)
            if isinstance(plot.get("spec"), dict):
                plot["spec"] = json.dumps(plot["spec"])
            spec["custom-plot"] = plot
        if self.detail:
            spec["display-mode"] = "detail"
        if self.hidden:
            spec["display-mode"] = "hidden"
        return spec


def ticks(label: str | None = None, detail: bool = False, precision: int = 3) -> Column:
    """A fraction, rate or measure as a tick on the range of the column's values, so the
    libraries compare at a glance; the number stays beside it."""
    plot = {"ticks": {"scale": "linear", "color": {"range": [INK]}}}
    return Column(label=label, precision=precision, plot=plot, detail=detail)


def bars(label: str | None = None, detail: bool = False, precision: int = 0) -> Column:
    """A count, with a bar proportional to its value."""
    plot = {"bars": {"scale": "linear", "color": {"range": [INK]}}}
    return Column(label=label, precision=precision, plot=plot, detail=detail)


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
    """A chart of the site from a long frame, one mark per row.

    ``line`` draws one thin line per ``series`` value over ``x``; ``bar`` draws bars of a
    nominal ``x`` side by side per series; ``dots`` draws one dot per row at ``x`` on the nominal
    ``y``; ``scatter`` draws one dot per row at the quantitative ``x`` and ``y``; ``strips``
    draws one tick per row at ``x``, in a panel per ``facet`` value with its own x scale (the
    strips of a summary's metrics).

    Few series (:data:`MAX_COLORED`) get a color each and a legend that toggles them; more share
    one ink, and a ``line`` or ``dots`` chart of many series becomes a heatmap (one row per
    series, the bins as columns, the value as color), which shows hundreds of libraries at once;
    ``many="ink"`` keeps the lines. Continuous axes zoom and pan; hovering a mark highlights its
    series, also across the panels of a strip chart; a mark links to the row of its series in
    the table view ``link``.
    """

    name: str
    title: str
    frame: pl.DataFrame
    x: str
    y: str | None = None
    desc: str = ""
    mark: str = "line"
    series: str = "library"
    #: ``quantitative`` or ``nominal`` (categories on the x axis)
    x_type: str = "quantitative"
    #: Column whose values become panels, side by side (``line``, ``bar``) or stacked
    #: (``strips``: one row per metric with its own x scale)
    facet: str | None = None
    #: The order of the panels; the order of the data without it
    facet_order: list[str] | None = None
    log_y: bool = False
    x_title: str | None = None
    y_title: str | None = None
    #: How many series are drawn: ``heatmap`` or ``ink`` (thin lines or dots in one ink)
    many: str = "heatmap"
    #: The table view whose rows the series link to (``None``: no links)
    link: str | None = "summary"
    #: Further columns that the tooltip of a mark names
    tooltip: tuple[str, ...] = ()

    @property
    def colored(self) -> bool:
        """Whether the series are few enough for a color each."""
        if self.series not in self.frame.columns:
            return False
        return self.frame[self.series].n_unique() <= MAX_COLORED

    @property
    def heatmap(self) -> bool:
        """Whether the chart is drawn as a heatmap: many series of a line or dots chart."""
        return self.mark in ("line", "dots") and not self.colored and self.many == "heatmap"

    def _color(self) -> dict:
        if self.colored:
            color: dict = {
                "field": self.series,
                "type": "nominal",
                "scale": {"range": list(COLORS)},
                "sort": None,
            }
            if self.frame[self.series].n_unique() == 1:  # the title names the one series
                color["legend"] = None
            return color
        return {"value": INK}

    def _tooltip(self) -> list[dict]:
        fields: list[dict] = [{"field": self.series, "type": "nominal"}]
        if self.facet:
            fields.append({"field": self.facet, "type": "nominal"})
        x: dict = {"field": self.x, "type": self.x_type}
        if self.x_type == "quantitative":
            x["format"] = ".4~g"
        fields.append(x)
        if self.y and self.mark != "dots":
            fields.append({"field": self.y, "type": "quantitative", "format": ".4~g"})
        elif self.y:
            fields.append({"field": self.y, "type": "nominal"})
        for extra in self.tooltip:
            fields.append({"field": extra, "type": "nominal"})
        return fields

    def _picked(self, base: float, hover: bool) -> tuple[dict, dict]:
        """The selection of a series (by the legend, or by hovering a mark) and the opacity that
        dims the other series while one is picked."""
        select: dict = {"type": "point", "fields": [self.series]}
        if hover:
            select |= {"on": "pointerover", "clear": "pointerout"}
        param = {"name": "pick", "select": select}
        if not hover:
            param["bind"] = "legend"
        opacity = {"condition": {"param": "pick", "value": base}, "value": round(base * 0.25, 2)}
        return param, opacity

    def _heatmap(self, encoding: dict, params: list) -> tuple[dict, int, int]:
        """The rect mark and axes of the heatmap: the series as rows, the x (or the category of a
        dots chart) as columns, the value as a sequential color."""
        if self.mark == "dots":
            column, column_title, value, value_title = self.y, self.y_title, self.x, self.x_title
            column_type = "nominal"
        else:
            column, column_title, value, value_title = self.x, self.x_title, self.y, self.y_title
            column_type = "ordinal" if self.x_type == "quantitative" else "nominal"
        axis = {"labelOverlap": "parity", "labelAngle": 0, "labelSeparation": 6}
        if column_type == "ordinal":
            axis["format"] = ".3~g"
        encoding["x"] = {
            "field": column,
            "type": column_type,
            "title": column if column_title is None else (column_title or None),
            "sort": None,
            "axis": axis,
        }
        encoding["y"] = {"field": self.series, "type": "nominal", "sort": None, "title": None}
        scale: dict = {"scheme": "blues", "domainMin": 0}  # a count or rate starts at zero
        if self.log_y:
            scale["type"] = "symlog"
        encoding["color"] = {
            "field": value,
            "type": "quantitative",
            "scale": scale,
            "legend": {"title": value_title or value, "gradientLength": 120},
        }
        encoding["tooltip"] = [
            {"field": self.series, "type": "nominal"},
            {"field": column, "type": column_type}
            | ({"format": ".4~g"} if column_type == "ordinal" else {}),
            {"field": value, "type": "quantitative", "format": ".4~g"},
        ]
        rows = int(self.frame[self.series].n_unique())
        return {"type": "rect", "tooltip": True}, 560, max(60, 12 * rows)

    def spec(self) -> dict:
        """The vega-lite specification of the chart."""
        x_title = self.x if self.x_title is None else (self.x_title or None)
        x: dict = {"field": self.x, "type": self.x_type, "title": x_title}
        encoding: dict = {"x": x, "color": self._color(), "tooltip": self._tooltip()}
        params: list[dict] = []
        mark: dict = {"type": self.mark, "tooltip": True}
        width, height = 560, 240
        zoom = {"name": "zoom", "select": "interval", "bind": "scales"}
        if self.heatmap:
            mark, width, height = self._heatmap(encoding, params)
        elif self.mark == "line":
            mark |= {"strokeWidth": 1.25 if self.colored else 1}
            encoding["detail"] = {"field": self.series, "type": "nominal"}
            param, encoding["opacity"] = self._picked(1 if self.colored else 0.45, not self.colored)
            params += [zoom, param]
        elif self.mark == "bar":
            mark |= {"cornerRadiusEnd": 2, "width": {"band": 0.6}}
            if self.x_type == "nominal":
                encoding["xOffset"] = {"field": self.series, "type": "nominal"}
            param, encoding["opacity"] = self._picked(0.9, not self.colored)
            params.append(param)
        elif self.mark == "dots":
            mark = {"type": "circle", "tooltip": True, "size": 36}
            encoding["y"] = {
                "field": self.y,
                "type": "nominal",
                "sort": None,
                "title": self.y_title or self.y,
            }
            height = max(120, 14 * int(self.frame[self.y].n_unique()))
            param, encoding["opacity"] = self._picked(
                0.8 if self.colored else 0.45, not self.colored
            )
            params += [zoom, param]
        elif self.mark == "scatter":
            mark = {"type": "circle", "tooltip": True, "size": 40}
            x["scale"] = {"padding": 12}  # the dots at the ends stay inside the axes
            encoding["y"] = {
                "field": self.y,
                "type": "quantitative",
                "title": self.y_title or self.y,
                "scale": {"padding": 12},
            }
            param, encoding["opacity"] = self._picked(
                0.8 if self.colored else 0.5, not self.colored
            )
            params += [zoom, param]
            width, height = 320, 320
        elif self.mark == "strips":
            mark = {"type": "tick", "tooltip": True, "thickness": 1.5}
            x["axis"] = {"format": ".3~g", "tickCount": 4}
            x["scale"] = {"zero": False, "nice": False, "padding": 8}
            width, height = 400, 28
            param, encoding["opacity"] = self._picked(0.8 if self.colored else 0.6, True)
            params.append(param)
        if self.mark in ("line", "bar") and not self.heatmap:
            y: dict = {"field": self.y, "type": "quantitative", "title": self.y_title or self.y}
            if self.log_y:
                y["scale"] = {"type": "log"}
            encoding["y"] = y
        if self.link and self.series in self.frame.columns:
            encoding["href"] = {"field": LINK}
        spec: dict = {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "title": self.title,
            "mark": mark,
            "encoding": encoding,
            "width": width,
            "height": height,
            "config": CHART_CONFIG,
        }
        if params:
            spec["params"] = params
        if self.facet:
            order = self.facet_order or self.frame[self.facet].unique(maintain_order=True).to_list()
            panel = {"field": self.facet, "type": "nominal", "sort": order, "title": None}
        if self.facet and self.mark == "strips":
            encoding["row"] = panel | {
                "header": {"labelAngle": 0, "labelAlign": "left", "labelAnchor": "start"}
            }
            spec["resolve"] = {"scale": {"x": "independent"}}
        elif self.facet:
            encoding["column"] = panel
            spec["width"] = 320
        return spec


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
    tables = {view.name: view for view in views if isinstance(view, Table)}
    for view in views:
        if not isinstance(view, Chart) or not view.link or view.link not in tables:
            continue
        if view.series in view.frame.columns and view.series in tables[view.link].frame.columns:
            datasets[view.name]["links"] = {
                LINK: {
                    "column": view.series,
                    "table-row": f"{view.link}/{view.series}",
                    "optional": True,
                }
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


def metric_strips(
    frame: pl.DataFrame, metrics: dict[str, str], name: str, title: str, desc: str = ""
) -> Chart:
    """A chart of the ``metrics`` of a summary ``frame`` (column -> label): one strip of ticks per
    metric, one tick per library, each strip on its own scale, so an outlying library shows."""
    columns = [c for c in metrics if c in frame.columns]
    long = (
        frame.select("library", *columns)
        .with_columns(pl.col(c).cast(pl.Float64, strict=False) for c in columns)
        .unpivot(index="library", on=columns, variable_name="metric", value_name="value")
        .with_columns(pl.col("metric").replace_strict(metrics, default=None))
        .drop_nulls("value")
    )
    return Chart(
        name,
        title,
        long,
        x="value",
        desc=desc,
        mark="strips",
        facet="metric",
        facet_order=[metrics[c] for c in columns],
        x_title="",
    )


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
    width = nice_width((high - low) / MAX_POINTS)
    return (
        frame.with_columns(((pl.col(x) / width).floor() * width).alias(x))
        .group_by([*groups, x], maintain_order=True)
        .agg(pl.col(y).sum())
    )


def nice_width(width: float) -> float:
    """The smallest of 1, 2 or 5 times a power of ten that is at least ``width``, so that bins
    start at round numbers."""
    if width <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(width))
    for factor in (1, 2, 5, 10):
        if factor * magnitude >= width:
            return float(factor * magnitude)
    return float(10 * magnitude)


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
