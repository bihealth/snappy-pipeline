# -*- coding: utf-8 -*-
"""Tests of the peddy summary: one row per sample with the sex and heterozygosity checks, the
pairs of every pedigree with their relatedness, and the scatter chart of the pairs."""

import json
from pathlib import Path

import fastjsonschema
import polars as pl

from snappy_wrappers.tools.qc_summary import Chart, config, summarize

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
KINDS = ("ped_check", "sex_check", "het_check")


def _trio():
    return summarize(
        "peddy",
        {kind: [str(FIXTURES / f"trio.{kind}.csv")] for kind in KINDS},
        ["P001-N1-DNA1-WES1"],
    )


def test_summary_joins_the_sex_and_heterozygosity_checks_per_sample():
    summary = _trio()
    frame = summary.table.frame

    assert frame.columns[:5] == ["library", "pedigree", "ped_sex", "predicted_sex", "sex_error"]
    assert frame["library"].to_list() == [
        "P001-N1-DNA1-WES1",
        "P002-N1-DNA1-WES1",
        "P003-N1-DNA1-WES1",
    ]
    assert frame["pedigree"].unique().to_list() == ["P001-N1-DNA1-WES1"]
    assert frame["predicted_sex"].to_list() == ["male", "male", "female"]
    assert frame["depth_outlier"].to_list() == ["False", "False", "True"]
    assert frame["mean_depth"].to_list() == [48.2, 51.0, 22.5]
    assert frame["sex_het_ratio"][2] == 1.107  # the chrX ratio of the sex check, kept apart
    assert frame["het_ratio"][2] == 0.426
    assert summary.table.columns["het_ratio"].plot["ticks"]
    assert summary.table.columns["sampled_sites"].detail


def test_pairs_table_and_scatter_compare_the_relatedness_of_the_ped_file_and_the_genotypes():
    views = {view.name: view for view in _trio().views}

    pairs = views["pairs"].frame
    assert pairs.columns[:6] == [
        "pedigree",
        "sample_a",
        "sample_b",
        "pedigree_relatedness",
        "rel",
        "rel_difference",
    ]
    assert pairs["rel"].to_list() == [0.49, 0.47, -0.02]
    assert pairs["pedigree_parents"].to_list() == ["True", "True", "False"]
    assert "pair" not in pairs.columns
    scatter = views["relatedness"]
    assert scatter.mark == "scatter" and scatter.series == "pedigree"
    spec = scatter.spec()
    assert spec["mark"]["type"] == "circle"
    assert spec["encoding"]["x"]["field"] == "pedigree_relatedness"
    y = spec["encoding"]["y"]
    assert (y["field"], y["type"], y["title"]) == (
        "rel",
        "quantitative",
        "relatedness from the genotypes",
    )
    assert y["scale"]["padding"] == 12  # the dots at the ends stay inside the axes
    assert {f["field"] for f in spec["encoding"]["tooltip"]} >= {"pedigree", "pair", "rel"}
    assert "href" not in spec["encoding"]  # pairs have no row in the summary table
    assert [p["name"] for p in spec["params"]] == ["zoom", "pick"]
    assert views["libraries"].frame["metric"].unique().sort().to_list() == [
        "IDR of the B allele frequency",
        "call rate",
        "chrX heterozygosity ratio",
        "heterozygosity ratio",
        "mean depth",
    ]


def test_the_configuration_of_the_peddy_site_is_valid():
    summary = _trio()
    cfg = config([summary.table, *summary.views], "checks: Pedigree checks", "/data")

    assert cfg["views"]["pairs"]["render-table"]["columns"]["pedigree"] == {
        "display-mode": "pinned"
    }
    assert "links" in cfg["datasets"]["libraries"] and "links" not in cfg["datasets"]["relatedness"]
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)


def test_summary_without_pedigrees_is_empty():
    summary = summarize("peddy", {kind: [] for kind in KINDS}, [])

    assert summary.table.frame.height == 0
    assert isinstance(summary.table.frame, pl.DataFrame)
    assert all(isinstance(view, (Chart, type(summary.table))) for view in summary.views)
