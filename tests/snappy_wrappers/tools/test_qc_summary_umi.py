# -*- coding: utf-8 -*-
"""Tests of the UMI consensus summary: one row per library with the samtools and fgumi columns,
the views of duplex and simplex consensus, the optional kinds and the datavzrd configuration."""

import json
from pathlib import Path

import fastjsonschema
import polars as pl
import pytest

from snappy_wrappers.tools.qc_summary import config, summarize, umi

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
LABELS = ["lib1", "lib2"]

#: fgumi files of two libraries by kind, the header first;
#: the test writes the fields of a row tab-separated
FILES = {
    "grouping_metrics": (
        [
            "accepted_sam_records discarded_non_pf discarded_poor_alignment "
            "discarded_ns_in_umi discarded_umis_to_short",
            "900 0 80 20 0",
        ],
        [
            "accepted_sam_records discarded_non_pf discarded_poor_alignment "
            "discarded_ns_in_umi discarded_umis_to_short",
            "950 0 50 0 0",
        ],
    ),
    "group_family_sizes": (
        [
            "family_size count fraction fraction_gt_or_eq_family_size",
            "1 60 0.6 1",
            "2 30 0.3 0.4",
            "40 10 0.1 0.1",
        ],
        [
            "family_size count fraction fraction_gt_or_eq_family_size",
            "1 50 0.5 1",
            "3 50 0.5 0.5",
        ],
    ),
    "position_group_sizes": (
        [
            "position_group_size count fraction fraction_gt_or_eq_position_group_size",
            "1 80 0.8 1",
            "2 20 0.2 0.2",
        ],
    )
    * 2,
    "umi_counts": (
        [
            "umi raw_observations raw_observations_with_errors unique_observations "
            "fraction_raw_observations fraction_unique_observations",
            "AAAA 400 0 200 0.4 0.4",
            "CCCC 600 10 300 0.6 0.6",
        ],
    )
    * 2,
    "consensus_stats": (
        [
            "key value description",
            "raw_reads_considered 1000 Total_raw_reads",
            "raw_reads_rejected 100 Rejected",
            "raw_reads_used 900 Used",
            "frac_raw_reads_used 0.9 Fraction_used",
            "raw_reads_rejected_for_insufficient_support 100 Insufficient_support",
            "consensus_reads_emitted 200 Emitted",
        ],
        [
            "key value description",
            "raw_reads_considered 500 Total_raw_reads",
            "raw_reads_rejected 0 Rejected",
            "raw_reads_used 500 Used",
            "frac_raw_reads_used 1 Fraction_used",
            "consensus_reads_emitted 100 Emitted",
        ],
    ),
    "filter_stats": (
        ["total_reads passed_reads failed_reads pass_rate", "200 180 20 0.9"],
        ["total_reads passed_reads failed_reads pass_rate", "100 100 0 1"],
    ),
    "correct_metrics": (
        [
            "umi total_matches perfect_matches one_mismatch_matches two_mismatch_matches "
            "other_matches fraction_of_matches representation",
            "AAAA 500 450 40 10 0 0.5 1.25",
            "CCCC 300 290 10 0 0 0.3 0.75",
            "NNNN 200 0 0 0 0 0.2 0.5",
        ],
    )
    * 2,
    "clip_metrics": (
        [
            "read_type reads reads_clipped_overlapping bases bases_clipped_overlapping",
            "Fragment 0 0 0 0",
            "ReadOne 100 10 9000 300",
            "ReadTwo 100 10 9000 300",
            "Pair 100 10 18000 600",
            "All 200 20 18000 600",
        ],
    )
    * 2,
}

#: The files of the consensus caller, by caller
CONSENSUS_FILES = {
    "duplex": {
        "family_sizes": (
            [
                "family_size cs_count cs_fraction cs_fraction_gt_or_eq_size ss_count ss_fraction "
                "ss_fraction_gt_or_eq_size ds_count ds_fraction ds_fraction_gt_or_eq_size",
                "1 50 0.5 1 60 0.6 1 40 0.4 1",
                "2 40 0.4 0.5 30 0.3 0.4 50 0.5 0.6",
                "31 10 0.1 0.1 10 0.1 0.1 10 0.1 0.1",
            ],
        )
        * 2,
        "duplex_family_sizes": (
            [
                "ab_size ba_size count fraction fraction_gt_or_eq_size",
                "1 0 50 0.5 1",
                "1 1 20 0.2 0.5",
                "2 1 20 0.2 0.3",
                "13 2 10 0.1 0.1",
            ],
        )
        * 2,
        "yield_metrics": (
            [
                "fraction read_pairs cs_families ss_families ds_families ds_duplexes "
                "ds_fraction_duplexes ds_fraction_duplexes_ideal",
                "0.5 500 400 450 420 50 0.119 0.2",
                "1 1000 700 800 750 150 0.2 0.4",
            ],
            [
                "fraction read_pairs cs_families ss_families ds_families ds_duplexes "
                "ds_fraction_duplexes ds_fraction_duplexes_ideal",
                "0.5 250 200 220 210 30 0.143 0.2",
                "1 500 350 400 375 75 0.2 0.25",
            ],
        ),
    },
    "simplex": {
        "family_sizes": (
            [
                "family_size cs_count cs_fraction cs_fraction_gt_or_eq_size ss_count ss_fraction "
                "ss_fraction_gt_or_eq_size",
                "1 50 0.5 1 60 0.6 1",
                "2 50 0.5 0.5 40 0.4 0.4",
            ],
        )
        * 2,
        "yield_metrics": (
            [
                "fraction read_pairs cs_families ss_families mean_ss_family_size ss_singletons "
                "ss_singleton_fraction ss_consensus_families",
                "0.5 500 400 450 1.1 300 0.667 150",
                "1 1000 700 800 1.25 400 0.5 400",
            ],
        )
        * 2,
    },
}


def _files(tmp_path, consensus, optional=("correct_metrics", "clip_metrics")):
    """The files of two libraries by kind: samtools from the fixtures, fgumi written to ``tmp_path``;
    ``optional`` names the optional kinds the task produces."""
    files = {
        kind: [str(FIXTURES / f"{label}.{kind}.txt") for label in LABELS]
        for kind in ("bamstats", "flagstats", "idxstats")
    }
    kinds = {k: v for k, v in FILES.items() if k not in ("correct_metrics", "clip_metrics")}
    kinds |= {k: FILES[k] for k in optional}
    for kind, contents in (kinds | CONSENSUS_FILES[consensus]).items():
        files[kind] = []
        for label, rows in zip(LABELS, contents):
            path = tmp_path / f"{label}.{kind}.txt"
            path.write_text("".join("\t".join(row.split()) + "\n" for row in rows))
            files[kind].append(str(path))
    return files


def _views(summary):
    return {view.name: view for view in summary.views}


def test_duplex_summary_has_one_row_per_library_with_the_umi_fractions(tmp_path):
    summary = summarize("umi", _files(tmp_path, "duplex"), LABELS, consensus="duplex")
    frame = summary.table.frame
    columns = summary.table.columns

    assert frame["library"].to_list() == LABELS
    assert frame.columns[:3] == ["library", "reads", "mapped_fraction"]  # samtools first
    assert frame["reads"].to_list() == [58722, 100000]
    assert columns["reads"].label == "consensus reads"
    assert frame["accepted_sam_records"].to_list() == [900, 950]
    assert frame["discarded_fraction"].to_list() == pytest.approx([0.1, 0.05])
    assert frame["mean_family_size"].to_list() == pytest.approx([5.2, 2.0])
    assert frame["frac_raw_reads_used"].to_list() == pytest.approx([0.9, 1.0])
    assert frame["consensus_reads_emitted"].to_list() == [200, 100]
    assert frame["pass_rate"].to_list() == pytest.approx([0.9, 1.0])
    # the yield of the full data
    assert frame["read_pairs"].to_list() == [1000, 500]
    assert frame["ds_duplexes"].to_list() == [150, 75]
    assert frame["ds_fraction_duplexes"].to_list() == pytest.approx([0.2, 0.2])
    assert frame["umi_corrected_fraction"].to_list() == pytest.approx([0.06, 0.06])
    assert frame["umi_unmatched_fraction"].to_list() == pytest.approx([0.2, 0.2])
    # a reason that one library lacks is empty in the other
    assert frame["raw_reads_rejected_for_insufficient_support"].to_list() == [100, None]
    assert columns["ds_duplexes"].plot["bars"] and columns["pass_rate"].plot["ticks"]
    assert columns["raw_reads_rejected_for_insufficient_support"].detail
    assert columns["discarded_poor_alignment"].detail and columns["ds_families"].detail
    assert "library" not in columns


def test_duplex_views_include_the_family_size_heatmap(tmp_path):
    views = _views(summarize("umi", _files(tmp_path, "duplex"), LABELS, consensus="duplex"))

    assert list(views) == [
        "libraries",
        "family_sizes",
        "group_family_sizes",
        "position_group_sizes",
        "duplex_family_sizes",
        "yield",
        "duplex_yield",
        "umi_counts",
        "umi_correction",
        "clipping",
        "insert_size",
        "coverage",
        "base_quality",
        "contigs",
        "idxstats",
    ]
    families = views["family_sizes"]
    assert families.facet == "kind" and families.facet_order == [
        "CS families",
        "SS families",
        "DS families",
    ]
    assert families.frame["family_size"].max() == 2  # sizes above 30 are left out
    assert views["group_family_sizes"].frame["family_size"].to_list() == [1, 2, 1, 3]
    heatmap = views["duplex_family_sizes"]
    assert heatmap.heatmap and heatmap.facet == "library" and heatmap.link is None
    lib1 = heatmap.frame.filter(pl.col("library") == "lib1")
    assert lib1.height == sum(ab + 1 for ab in range(1, 13))  # the triangle up to AB 12
    assert lib1["ab"].unique(maintain_order=True).to_list()[:2] == ["AB 1", "AB 2"]
    assert lib1.filter(pl.col("ab") == "AB 2")["families"].to_list() == [0, 20, 0]
    assert lib1["families"].sum() == 90  # AB 13 is left out
    spec = heatmap.spec()
    assert spec["mark"]["type"] == "rect" and spec["encoding"]["column"]["field"] == "library"
    assert spec["encoding"]["color"]["scale"]["type"] == "symlog"
    assert views["yield"].y == "ds_duplexes"
    ratio = views["duplex_yield"].frame
    assert ratio.filter(pl.col("fraction") == 1)["of_ideal"].to_list() == pytest.approx([0.5, 0.8])
    metrics = views["libraries"].frame["metric"].unique().to_list()
    assert {"duplex fraction", "corrected UMIs", "mapped"} <= set(metrics)
    assert views["umi_counts"].frame["umi"].to_list() == ["CCCC", "AAAA"] * 2  # most observed first
    assert views["clipping"].frame.height == 10


def test_umi_counts_keep_the_most_observed_umis(tmp_path, monkeypatch):
    monkeypatch.setattr(umi, "MAX_UMIS", 1)

    summary = summarize("umi", _files(tmp_path, "duplex"), LABELS, consensus="duplex")

    assert _views(summary)["umi_counts"].frame["umi"].to_list() == ["CCCC", "CCCC"]


def test_simplex_summary_without_the_optional_kinds(tmp_path):
    files = _files(tmp_path, "simplex", optional=())
    assert "duplex_family_sizes" not in files

    summary = summarize("umi", files, LABELS, consensus="simplex")

    frame = summary.table.frame
    assert frame["ss_consensus_families"].to_list() == [400, 400]
    assert frame["ss_singleton_fraction"].to_list() == pytest.approx([0.5, 0.5])
    assert "ds_duplexes" not in frame.columns
    assert "umi_corrected_fraction" not in frame.columns
    views = _views(summary)
    assert not {"duplex_family_sizes", "duplex_yield", "umi_correction", "clipping"} & set(views)
    assert views["family_sizes"].frame["kind"].unique().sort().to_list() == [
        "CS families",
        "SS families",
    ]
    assert views["yield"].y == "ss_consensus_families"
    assert "SS singletons" in views["libraries"].frame["metric"].to_list()


@pytest.mark.parametrize("consensus", ["duplex", "simplex"])
def test_the_configuration_of_the_umi_site_is_valid(tmp_path, consensus):
    summary = summarize("umi", _files(tmp_path, consensus), LABELS, consensus=consensus)

    cfg = config([summary.table, *summary.views], "umi: UMI QC", "/data")

    assert cfg["default-view"] == "summary"
    assert cfg["datasets"]["family_sizes"]["links"]["row"]["table-row"] == "summary/library"
    columns = cfg["views"]["umi_counts"]["render-table"]["columns"]
    assert columns["umi"] == {"display-mode": "pinned"}
    if consensus == "duplex":
        assert "links" not in cfg["datasets"]["duplex_family_sizes"]
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)  # fills in the schema's defaults
