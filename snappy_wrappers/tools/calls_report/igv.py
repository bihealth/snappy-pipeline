# -*- coding: utf-8 -*-
"""IGV views of the top calls with igv-reports.

The top calls of the chosen views, taken round-robin by rank so that every view contributes,
become the sites of one ``create_report`` run:
a page with an igv.js session per site, the reads of every library sliced around it,
and a table of the sites.
A structural variant spans its ``END``; one longer than ``max_span`` gets two sites,
its start and its end breakpoint, so one large deletion cannot fill the page.
The page is cut to a size limit by dropping the lowest-ranked sites and rerunning.
igv.js is inlined from the copy that the environment's post-deploy script fetched,
so no job needs network,
and the calls table gets an ``igv`` column with the rank of a selected call, linked to the page.
The link carries the rank in the URL hash, which a script in the page resolves to the call's row;
the template of igv-reports preselects the first row only.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import Column, Summary, Table
from snappy_wrappers.tools.vcf_events import select_names

#: The page inside the site, relative to its root
PAGE = "igv/top.html"

#: The script tag of igv-reports' template, replaced by the inlined igv.js
IGV_JS_TAG = re.compile(
    r'<script[^>]*src="https://cdn\.jsdelivr\.net/npm/igv@3\.5\.2/dist/igv\.min\.js"[^>]*>\s*</script>'
)

#: Columns of the sites file, in order; ``create_report`` reads the sequence, begin and end;
#: ``breakpoint`` is ``start`` or ``end`` for the two sites of a long structural variant
SITE_COLUMNS = (
    "rank",
    "view",
    "chromosome",
    "begin",
    "end",
    "breakpoint",
    "ref",
    "alt",
    "posterior",
    "gene",
    "consequence",
)

KEY = ("chromosome", "position", "ref", "alt")

#: Selects the site of the URL hash, the rank of a call as the ``igv`` column links it;
#: the template builds its table once igv.js is up, so the script waits for the table
HASH_SCRIPT = """<script type="text/javascript">
    (function () {
        function select() {
            const table = document.getElementById("variant_table")
            if (!table) {
                window.setTimeout(select, 100)
                return
            }
            const rank = decodeURIComponent(window.location.hash.slice(1))
            if (!rank) return
            const headers = Array.from(table.querySelectorAll("th"), (th) => th.textContent)
            const column = headers.indexOf("rank")
            for (const row of table.querySelectorAll("tbody tr")) {
                const cell = row.cells[column]
                if (cell && cell.textContent === rank) {
                    cell.click()  // the template's handler loads the session of the row
                    row.scrollIntoView({block: "center"})
                    return
                }
            }
        }
        window.addEventListener("hashchange", select)
        select()
    })()
</script>
"""


def igv_js_path() -> Path:
    """The igv.js that the environment's post-deploy script fetched."""
    return Path(sys.prefix) / "share" / "igv-reports" / "igv.min.js"


def chosen_views(views: Sequence[Table], names: Sequence[str]) -> list[Table]:
    """The views that ``names`` select:
    names, which must exist, or patterns such as ``denovo_*`` (:func:`select_names`);
    all views without names."""
    if not names:
        return list(views)
    selected = select_names(list(names), [view.name for view in views], what="view")
    return [view for view in views if view.name in selected]


def _site(
    rank: int, view: str, call: dict[str, Any], begin: int, end: int, breakpoint: str
) -> dict[str, Any]:
    return {
        "rank": rank,
        "view": view,
        "chromosome": call["chromosome"],
        "begin": begin,
        "end": end,
        "breakpoint": breakpoint,
        "ref": call["ref"],
        "alt": call["alt"],
        "posterior": call.get(f"P({view})"),
        "gene": call.get("gene", ""),
        "consequence": call.get("consequence", ""),
    }


def select_sites(
    views: Sequence[Table], top: int, names: Sequence[str] = (), *, max_span: int
) -> pl.DataFrame:
    """The sites of the page:
    the first ``top`` calls of each view in ``names`` (all views if empty),
    interleaved by rank so that every view contributes, each call once.
    A call spans to its ``end`` (a structural variant's ``END``) when the table has one;
    a span longer than ``max_span`` becomes two sites, the start and the end breakpoint."""
    chosen = chosen_views(views, names)
    heads = [view.frame.head(top).to_dicts() for view in chosen]
    rows: list[dict[str, Any]] = []
    seen: set[tuple] = set()
    for index in range(top):
        for view, head in zip(chosen, heads):
            if index >= len(head):
                continue
            call = head[index]
            key = tuple(call[column] for column in KEY)
            if key in seen:
                continue
            seen.add(key)
            begin = call["position"]
            end = call.get("end") or begin + len(call["ref"]) - 1
            if end - begin + 1 > max_span:
                rows.append(_site(len(seen), view.name, call, begin, begin, "start"))
                rows.append(_site(len(seen), view.name, call, end, end, "end"))
            else:
                rows.append(_site(len(seen), view.name, call, begin, end, ""))
    schema = {
        "rank": pl.Int64,
        "view": pl.String,
        "chromosome": pl.String,
        "begin": pl.Int64,
        "end": pl.Int64,
        "breakpoint": pl.String,
        "ref": pl.String,
        "alt": pl.String,
        "posterior": pl.Float64,
        "gene": pl.String,
        "consequence": pl.String,
    }
    return pl.DataFrame(rows, schema=schema)


def track_config(
    labels: Sequence[str], bams: Sequence[str], bais: Sequence[str], vcf: str, tbi: str
) -> list[dict[str, str]]:
    """The tracks of every session: the calls and one alignment track per library."""
    tracks = [
        {
            "name": "calls",
            "type": "variant",
            "format": "vcf",
            "url": os.path.abspath(vcf),
            "indexURL": os.path.abspath(tbi),
        }
    ]
    for label, bam, bai in zip(labels, bams, bais):
        tracks.append(
            {
                "name": label,
                "type": "alignment",
                "format": "bam",
                "url": os.path.abspath(bam),
                "indexURL": os.path.abspath(bai),
            }
        )
    return tracks


def create_report_command(
    sites: str,
    fasta: str,
    tracks: str,
    flanking: int,
    output: str,
    title: str,
    extra_arguments: Sequence[str] = (),
) -> list[str]:
    """``create_report`` on the sites file (1-based begin and end in its columns 4 and 5)."""
    return [
        "create_report",
        sites,
        "--fasta",
        fasta,
        "--track-config",
        tracks,
        "--sequence",
        str(SITE_COLUMNS.index("chromosome") + 1),
        "--begin",
        str(SITE_COLUMNS.index("begin") + 1),
        "--end",
        str(SITE_COLUMNS.index("end") + 1),
        "--flanking",
        str(flanking),
        "--title",
        title,
        "--output",
        output,
        *extra_arguments,
    ]


def fit_to_size(sites: int, size: int, limit: int) -> int:
    """How many of ``sites`` keep a page of ``size`` bytes under ``limit``:
    all if it fits, else the proportion with a tenth of margin, so that the next page fits."""
    if size <= limit:
        return sites
    return max(0, math.floor(sites * 0.9 * limit / size))


def inline_igv_js(html: str, js: str) -> str:
    """``html`` with the template's igv.js script tag replaced by the script itself."""
    if len(IGV_JS_TAG.findall(html)) != 1:
        raise ValueError(
            "the igv-reports page has no igv.js script tag for igv@3.5.2 to inline; the "
            "igv-reports version changed, update the post-deploy script"
        )
    return IGV_JS_TAG.sub(lambda _: f'<script type="text/javascript">{js}</script>', html, count=1)


def select_by_hash(html: str) -> str:
    """``html`` with :data:`HASH_SCRIPT` before its last ``</body>``
    (the inlined igv.js holds the tag in a string too)."""
    at = html.rfind("</body>")
    if at < 0:
        raise ValueError(
            "the igv-reports page has no </body> to add the hash selection to; the igv-reports "
            "version changed"
        )
    return html[:at] + HASH_SCRIPT + html[at:]


def build_page(
    summary: Summary,
    settings: dict[str, Any],
    bams: Sequence[str],
    bais: Sequence[str],
    vcf: str,
    tbi: str,
    fasta: str,
    workdir: str,
    title: str,
    run: Callable[..., Any] = subprocess.run,
) -> tuple[str | None, pl.DataFrame, int]:
    """Create the page in ``workdir``:
    the path (``None`` without sites),
    the sites it shows and the number of calls requested before the size limit cut them."""
    sites = select_sites(
        summary.views, settings["top"], settings.get("views", ()), max_span=settings["max_span"]
    )
    requested = sites["rank"].n_unique()
    limit = settings["max_size_mb"] * 1024 * 1024
    tracks_path = os.path.join(workdir, "tracks.json")
    with open(tracks_path, "wt") as handle:
        json.dump(track_config(settings["libraries"], bams, bais, vcf, tbi), handle)
    page = os.path.join(workdir, "top.html")
    while sites.height:
        sites_path = os.path.join(workdir, "sites.tsv")
        # Unquoted: igv-reports shows a quoted empty string as "" in its table
        sites.write_csv(sites_path, separator="\t", quote_style="never")
        command = create_report_command(
            sites_path,
            fasta,
            tracks_path,
            settings["flanking"],
            page,
            title,
            settings.get("extra_arguments", ()),
        )
        run(command, check=True)
        size = os.path.getsize(page)
        keep = fit_to_size(sites.height, size, limit)
        if keep == sites.height:
            break
        sites = sites.head(keep)
    if not sites.height:
        return None, sites, requested
    html = inline_igv_js(
        Path(page).read_text(encoding="utf-8"), igv_js_path().read_text(encoding="utf-8")
    )
    Path(page).write_text(select_by_hash(html), encoding="utf-8")
    return page, sites, requested


def annotate(summary: Summary, sites: pl.DataFrame, requested: int, limit_mb: int) -> None:
    """Add the ``igv`` column (the rank of a selected call,
    linked to the page and to IGV desktop)
    to every table of ``summary`` and a sentence to the descriptions."""
    calls = sites.filter(pl.col("breakpoint") != "end")  # one row per call
    ranks = calls.select(
        *KEY[:1], pl.col("begin").alias("position"), "ref", "alt", pl.col("rank").alias("igv")
    )
    # IGV desktop shows the whole structural variant
    locus = "{chromosome}:{position}"
    if "end" in summary.table.frame.columns:
        locus += "-{end}"
    column = Column(
        precision=0,
        links={
            "IGV": f"../{PAGE}#{{value}}",
            "IGV desktop": f"http://localhost:60151/goto?locus={locus}",
        },
    )
    for table in (summary.table, *summary.views):
        table.frame = table.frame.join(ranks, on=list(KEY), how="left", maintain_order="left")
        table.columns = {**table.columns, "igv": column}
    per_view = calls.group_by("view").len().to_dicts() if calls.height else []
    counts = {row["view"]: row["len"] for row in per_view}
    for view in summary.views:
        if view.name in counts:
            n = counts[view.name]
            view.desc += (
                f"\n- The top {n} {'call' if n == 1 else 'calls'} of this view "
                f"{'has' if n == 1 else 'have'} IGV views in `{PAGE}` (the `igv` column)."
            )
    note = (
        f"\n- IGV views of the top calls of the views are in `{PAGE}`; "
        "the `igv` column links a call to its view."
    )
    if requested > calls.height:
        note += (
            f"\n- {calls.height} of the {requested} requested calls fit the size limit of "
            f"{limit_mb} MB."
        )
    summary.table.desc += note


def install_page(page: str | None, site: str) -> None:
    """Move the page into the rendered site (after datavzrd wrote it)."""
    if page is None:
        return
    target = os.path.join(site, PAGE)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    shutil.move(page, target)
