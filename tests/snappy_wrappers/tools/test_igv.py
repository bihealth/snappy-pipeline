# -*- coding: utf-8 -*-
"""Tests of the IGV views of the calls report: the sites, the tracks, the command, the size
fit, the inlining of igv.js and the igv column. ``create_report`` itself is not run."""

import json
import os
from pathlib import Path

import polars as pl
import pytest

from snappy_wrappers.tools.calls_report import Options, igv, summarize
from snappy_wrappers.tools.qc_summary import Summary, Table

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"


def _view(name, positions, posteriors=None):
    frame = pl.DataFrame(
        {
            "chromosome": ["1"] * len(positions),
            "position": positions,
            "ref": ["A"] * len(positions),
            "alt": ["T"] * len(positions),
            f"P({name})": posteriors or [0.9] * len(positions),
            "gene": ["G"] * len(positions),
        }
    )
    return Table(name, name, frame)


def test_select_sites_interleaves_the_views_by_rank_and_takes_each_call_once():
    views = [_view("denovo", [10, 20, 30]), _view("present", [20, 40]), _view("dominant", [])]

    sites = igv.select_sites(views, top=2)

    # Rank 1 of every view, then rank 2; denovo's second call is present's first, counted once
    assert sites["rank"].to_list() == [1, 2, 3]
    assert list(zip(sites["view"].to_list(), sites["begin"].to_list())) == [
        ("denovo", 10),
        ("present", 20),
        ("present", 40),
    ]
    assert sites["end"].to_list() == [10, 20, 40]
    assert sites["posterior"].to_list() == [0.9, 0.9, 0.9]
    assert list(sites.columns) == list(igv.SITE_COLUMNS)
    # Only the named views, top above their rows takes all they have
    named = igv.select_sites(views, top=10, names=["present"])
    assert named["begin"].to_list() == [20, 40]
    assert igv.select_sites(views, top=5, names=["dominant"]).height == 0


def test_track_config_lists_the_calls_and_one_alignment_track_per_library():
    tracks = igv.track_config(
        ["lib1", "lib2"], ["a.bam", "b.bam"], ["a.bai", "b.bai"], "c.vcf.gz", "c.tbi"
    )

    assert [t["name"] for t in tracks] == ["calls", "lib1", "lib2"]
    assert tracks[0]["type"] == "variant" and tracks[0]["url"] == os.path.abspath("c.vcf.gz")
    assert tracks[1] == {
        "name": "lib1",
        "type": "alignment",
        "format": "bam",
        "url": os.path.abspath("a.bam"),
        "indexURL": os.path.abspath("a.bai"),
    }


def test_create_report_command_names_the_site_columns_one_based():
    command = igv.create_report_command(
        "sites.tsv", "ref.fa", "tracks.json", 300, "top.html", "fdr: calls", ["--subsample", "0.5"]
    )

    assert command[:2] == ["create_report", "sites.tsv"]
    assert command[command.index("--sequence") + 1] == "3"
    assert command[command.index("--begin") + 1] == "4"
    assert command[command.index("--end") + 1] == "5"
    assert command[command.index("--flanking") + 1] == "300"
    assert command[-2:] == ["--subsample", "0.5"]
    assert "--standalone" not in command


def test_fit_to_size_keeps_all_sites_that_fit_and_trims_with_a_margin():
    assert igv.fit_to_size(20, 10_000, 50_000) == 20
    assert igv.fit_to_size(20, 100_000, 50_000) == 9  # floor(20 * 0.9 * 0.5)
    assert igv.fit_to_size(1, 100_000, 50_000) == 0


def test_inline_igv_js_replaces_the_template_tag_once():
    html = '<head>\n<script src="https://cdn.jsdelivr.net/npm/igv@3.5.2/dist/igv.min.js"></script>\n</head>'

    inlined = igv.inline_igv_js(html, "var igv = 1;")

    assert inlined == '<head>\n<script type="text/javascript">var igv = 1;</script>\n</head>'
    with pytest.raises(ValueError, match="igv-reports version changed"):
        igv.inline_igv_js("<head></head>", "x")


def test_build_page_trims_the_sites_until_the_page_fits(tmp_path, monkeypatch):
    summary = summarize(str(FIXTURES / "calls.vcf"), Options())
    sizes = iter([30 * 1024 * 1024, 2 * 1024 * 1024])  # the first page is too big
    runs = []

    def run(command, check):
        runs.append(command)
        sites = pl.read_csv(command[1], separator="\t")
        Path(command[command.index("--output") + 1]).write_text(
            '<script src="https://cdn.jsdelivr.net/npm/igv@3.5.2/dist/igv.min.js"></script>'
            + "x" * next(sizes)
        )
        assert sites.height in (6, 1)

    monkeypatch.setattr(igv, "igv_js_path", lambda: tmp_path / "igv.min.js")
    (tmp_path / "igv.min.js").write_text("var igv;")
    settings = {
        "top": 10,
        "views": [],
        "flanking": 300,
        "max_size_mb": 10,
        "libraries": ["a", "b", "c"],
    }

    page, sites, requested = igv.build_page(
        summary,
        settings,
        ["a.bam", "b.bam", "c.bam"],
        ["a.bai", "b.bai", "c.bai"],
        "c.vcf.gz",
        "c.tbi",
        "ref.fa",
        str(tmp_path),
        "T",
        run=run,
    )

    assert requested == 6 and sites.height == 1  # floor(6 * 0.9 * 10 / 30)
    assert len(runs) == 2
    assert Path(page).read_text().startswith('<script type="text/javascript">var igv;</script>')
    tracks = json.loads((tmp_path / "tracks.json").read_text())
    assert [t["name"] for t in tracks] == ["calls", "a", "b", "c"]

    igv.annotate(summary, sites, requested, 10)

    assert summary.table.frame["igv"].to_list() == [1, None, None, None, None, None]
    present = next(view for view in summary.views if view.name == "present")
    assert present.frame["igv"].to_list()[0] == 1
    assert present.columns["igv"].links == {
        "IGV": "../igv/top.html#{value}",
        "IGV desktop": "http://localhost:60151/goto?locus={chromosome}:{position}",
    }
    assert "The top 1 call of this view has IGV views in igv/top.html" in present.desc
    assert "1 of the 6 requested calls fit the size limit of 10 MB" in summary.table.desc

    site = tmp_path / "site"
    site.mkdir()
    igv.install_page(page, str(site))
    assert (site / "igv" / "top.html").is_file()
    igv.install_page(None, str(site))  # nothing to install


def test_build_page_without_sites_writes_no_page(tmp_path):
    summary = Summary(table=_view("all", []), views=[_view("present", [])])
    settings = {"top": 5, "views": [], "flanking": 300, "max_size_mb": 50, "libraries": []}

    page, sites, requested = igv.build_page(
        summary, settings, [], [], "c.vcf.gz", "c.tbi", "ref.fa", str(tmp_path), "T", run=None
    )

    assert page is None and sites.height == 0 and requested == 0
