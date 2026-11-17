# -*- coding: utf-8 -*-
"""Tests of the QC summary tools: the parsers of each format, the datavzrd configuration and the
thinning of distributions."""

import json
import shutil
from pathlib import Path

import fastjsonschema
import polars as pl
import pytest

from snappy_wrappers.tools import qc_summary
from snappy_wrappers.tools.qc_summary import Chart, Table, config, numeric, summarize, thin

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"


def _validate_config(cfg):
    """``cfg`` is valid against datavzrd's configuration schema."""
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)


def _views(summary):
    return {view.name: view for view in summary.views}


def _validate(summary, title="QC"):
    """The datavzrd configuration of ``summary`` is valid against the schema."""
    cfg = config([summary.table, *summary.views], title, "/tmp")
    _validate_config(cfg)
    return cfg


def test_samtools_summary_has_one_row_per_library_with_the_fractions():
    summary = summarize(
        "samtools",
        {
            kind: [str(FIXTURES / f"lib{i}.{kind}.txt") for i in (1, 2)]
            for kind in ("bamstats", "flagstats", "idxstats")
        },
        ["lib1", "lib2"],
    )

    frame = summary.table.frame
    assert frame.columns[:3] == ["library", "reads", "mapped_fraction"]
    assert frame["library"].to_list() == ["lib1", "lib2"]
    assert frame["reads"].to_list() == [58722, 100000]
    assert frame["mapped_fraction"][0] == pytest.approx(58448 / 58722)
    assert frame["duplicated_fraction"].to_list() == pytest.approx([0.0, 0.05])
    assert frame["properly_paired_fraction"][0] == pytest.approx(57956 / 58722)
    assert frame["singletons"].to_list() == [274, 274]
    views = _views(summary)
    assert set(views) == {
        "libraries",
        "insert_size",
        "coverage",
        "base_quality",
        "contigs",
        "idxstats",
    }
    strips = views["libraries"]
    assert strips.mark == "strips" and strips.facet == "metric"
    assert strips.frame.filter(pl.col("metric") == "mapped")["value"].to_list() == pytest.approx(
        [58448 / 58722, 0.9]
    )
    assert views["contigs"].mark == "dots" and views["contigs"].y == "contig"
    assert views["base_quality"].frame["read"].unique().sort().to_list() == ["first", "last"]
    assert views["base_quality"].frame.height == 2 * 2 * 6
    assert views["insert_size"].frame.columns == ["library", "insert_size", "pairs"]
    idxstats = views["idxstats"].frame
    assert idxstats["contig"].to_list() == ["12", "12"]
    assert idxstats["mapped_per_mb"][1] == pytest.approx(90000 / 27500000 * 1e6)
    _validate(summary)


def test_alfred_summary_lists_the_fraction_of_targets_above_the_coverage_levels():
    summary = summarize(
        "alfred", {"json": [str(FIXTURES / "lib1.target_coverage.json.gz")]}, ["lib1"]
    )

    frame = summary.table.frame
    assert frame.columns[:4] == ["library", "read_group", "MappedFraction", "MappedProperFraction"]
    assert frame["library"].to_list() == ["lib1"]
    assert "#Mapped" not in frame.columns and frame["Mapped"][0] == 58448
    assert frame["OnTargetFraction"][0] == pytest.approx(3702647 / 5877555)
    assert frame["TargetsAbove20x"][0] == pytest.approx(0.058347, abs=1e-5)
    assert summary.table.columns["AlignedBases"].detail
    views = _views(summary)
    assert set(views) == {
        "libraries",
        "target_coverage",
        "on_target",
        "coverage",
        "insert_size",
        "base_quality",
        "gc_content",
        "contigs",
        "chromosomes",
    }
    assert "targets above 20x" in views["libraries"].frame["metric"].to_list()
    contigs = views["contigs"].frame
    assert contigs["mapped_per_mb"][0] == pytest.approx(
        contigs["mapped"][0] / (contigs["Size"][0] - contigs["N"][0]) * 1e6
    )
    assert views["base_quality"].frame["read"].unique().sort().to_list() == ["Read1", "Read2"]
    assert views["gc_content"].frame["library"].unique().sort().to_list() == [
        "lib1",
        "reference",
        "target",
    ]
    assert views["insert_size"].frame.height == 12  # the orientations summed
    assert views["chromosomes"].frame.columns[:3] == ["library", "Chr", "Size"]
    _validate(summary)


def test_picard_summary_joins_the_programs_into_one_row_per_library():
    programs = {
        "CollectHsMetrics": "CollectHsMetrics",
        "CollectAlignmentSummaryMetrics": "CollectMultipleMetrics.alignment_summary_metrics",
        "CollectInsertSizeMetrics": "CollectMultipleMetrics.insert_size_metrics",
    }
    summary = summarize(
        "picard",
        {
            program: [str(FIXTURES / f"lib{i}.{suffix}.txt") for i in (1, 2)]
            for program, suffix in programs.items()
        },
        ["lib1", "lib2"],
    )

    frame = summary.table.frame
    assert frame["library"].to_list() == ["lib1", "lib2"]
    assert frame["HsMetrics.MEAN_TARGET_COVERAGE"].to_list() == [80.5, 120.25]
    assert frame["AlignmentSummaryMetrics.CATEGORY"].to_list() == ["PAIR", "PAIR"]
    assert frame["AlignmentSummaryMetrics.TOTAL_READS"].to_list() == [100000, 200000]
    assert frame["InsertSizeMetrics.MEDIAN_INSERT_SIZE"].to_list() == [210, 230]
    assert summary.table.columns["HsMetrics.PCT_TARGET_BASES_20X"].plot["ticks"]
    views = _views(summary)
    assert set(views) == {
        "libraries",
        "HsMetrics",
        "AlignmentSummaryMetrics",
        "InsertSizeMetrics",
        "target_coverage",
        "insert_size",
    }
    metrics = views["libraries"].frame["metric"].unique().to_list()
    assert "pct target bases 20x" in metrics and "mean target coverage" in metrics
    assert views["AlignmentSummaryMetrics"].frame.height == 6  # three categories per library
    insert_size = views["insert_size"].frame
    assert insert_size.columns == ["library", "insert_size", "read pairs"]
    assert insert_size.filter(pl.col("library") == "lib1")["read pairs"].to_list() == [11, 52, 21]
    assert views["target_coverage"].log_y
    _validate(summary)


def test_bcftools_stats_summary_has_one_row_per_sample_and_variant_set():
    summary = summarize(
        "bcftools_stats", {"txt": [str(FIXTURES / "lib1.bcftools_stats.txt")]}, ["lib1"]
    )

    frame = summary.table.frame
    assert frame.columns[:4] == ["library", "set", "het", "hom_alt"]
    assert frame["set"].unique().sort().to_list() == ["known", "novel"]
    # The fixture is a trio VCF with three samples, one row each per set
    assert frame.height == 6
    novel = frame.filter(pl.col("set") == "novel")
    assert novel["het"].to_list() == [125, 125, 125]
    assert novel["ts_tv"][0] == pytest.approx(257 / 128)
    assert novel["het_hom_ratio"][0] == pytest.approx(125 / 260)
    assert novel["sites"][0] == 433
    views = _views(summary)
    assert set(views) == {"libraries", "substitutions", "indel_length", "depth", "quality"}
    assert views["substitutions"].mark == "bar" and views["substitutions"].x_type == "nominal"
    assert all(view.facet == "set" for name, view in views.items() if name != "libraries")
    # one tick per sample and variant set
    assert views["libraries"].frame["library"][0].endswith(" (known)")
    assert views["depth"].frame["bin"].dtype == pl.Float64
    _validate(summary)


def test_bcftools_stats_without_split_by_id_has_one_set(tmp_path):
    text = (FIXTURES / "lib1.bcftools_stats.txt").read_text(encoding="utf-8")
    lines = [line for line in text.splitlines() if "\t1\t" not in line[:8]]
    lines = [
        line.replace('calls.vcf.gz:known (sites with ID different from ".")', "calls.vcf.gz")
        for line in lines
    ]
    path = tmp_path / "single.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary = summarize("bcftools_stats", {"txt": [str(path)]}, ["lib1"])

    assert summary.table.frame["set"].unique().to_list() == ["all"]
    assert all(view.facet is None for view in summary.views if view.name != "libraries")


def test_thin_bins_long_distributions_per_library():
    frame = pl.DataFrame(
        {
            "library": ["a"] * 1000 + ["b"] * 10,
            "x": list(range(1000)) + list(range(10)),
            "y": [1] * 1010,
        }
    )

    thinned = thin(frame, "x", "y")

    assert thinned.filter(pl.col("library") == "a").height <= qc_summary.MAX_POINTS + 1
    assert thinned.filter(pl.col("library") == "a")["y"].sum() == 1000
    assert thinned.filter(pl.col("library") == "b")["y"].sum() == 10
    assert thin(frame.head(10), "x", "y").height == 10  # short distributions stay


def test_numeric_casts_only_columns_of_numbers():
    frame = pl.DataFrame({"n": ["1", "2"], "f": ["1.5", ""], "t": ["a", "1"]})

    result = numeric(frame)

    assert result["n"].dtype == pl.Int64 and result["n"].to_list() == [1, 2]
    assert result["f"].dtype == pl.Float64 and result["f"].to_list() == [1.5, None]
    assert result["t"].dtype == pl.String


def test_config_pins_the_first_columns_and_plots_the_charts():
    table = Table(
        "summary",
        "T",
        pl.DataFrame({"library": ["a"], "v": [0.5]}),
        columns={"v": qc_summary.ticks()},
    )
    chart = Chart(
        "c",
        "C",
        pl.DataFrame({"library": ["a", "a"], "read": ["first", "last"], "x": [1, 1], "y": [2, 3]}),
        x="x",
        y="y",
        facet="read",
    )

    cfg = config([table, chart], "task: T", "/data")

    assert cfg["default-view"] == "summary"
    assert cfg["datasets"]["c"] == {"path": "/data/c.tsv", "separator": "\t"}
    columns = cfg["views"]["summary"]["render-table"]["columns"]
    assert columns["library"] == {"display-mode": "pinned"}
    assert columns["v"]["plot"]["ticks"]["scale"] == "linear"
    spec = json.loads(cfg["views"]["c"]["render-plot"]["spec"])
    assert spec["encoding"]["column"]["field"] == "read"
    assert spec["encoding"]["column"]["sort"] == ["first", "last"]
    assert spec["encoding"]["color"]["field"] == "library"
    assert spec["encoding"]["color"]["legend"] is None  # one series: the title names it
    assert spec["config"]["axis"]["grid"] is False and spec["config"]["view"]["stroke"] is None
    _validate_config(cfg)


def test_many_libraries_share_one_ink_and_few_get_colors():
    few = pl.DataFrame({"library": ["a", "b"], "x": [1, 2], "y": [3, 4]})
    many = pl.DataFrame({"library": [f"l{i}" for i in range(7)], "x": range(7), "y": range(7)})

    colored = Chart("c", "C", few, x="x", y="y").spec()
    plain = Chart("c", "C", many, x="x", y="y").spec()

    assert colored["encoding"]["color"]["field"] == "library"
    assert colored["encoding"]["color"]["scale"]["range"][:2] == list(qc_summary.COLORS[:2])
    assert plain["encoding"]["color"] == {"value": qc_summary.INK}
    assert plain["encoding"]["detail"]["field"] == "library"
    assert plain["mark"]["opacity"] < 1 and plain["mark"]["strokeWidth"] <= 2


def test_metric_strips_put_each_metric_on_its_own_scale():
    frame = pl.DataFrame({"library": ["a", "b"], "f": [0.5, 0.6], "n": [10, 20], "t": ["x", "y"]})

    chart = qc_summary.metric_strips(frame, {"f": "fraction", "n": "count", "t": "text"}, "m", "M")
    spec = chart.spec()

    assert chart.frame["metric"].to_list() == ["fraction", "fraction", "count", "count"]
    assert spec["mark"]["type"] == "tick"
    assert spec["encoding"]["row"]["field"] == "metric"
    assert spec["resolve"] == {"scale": {"x": "independent"}}
    assert spec["encoding"]["x"]["title"] is None


@pytest.mark.skipif(shutil.which("datavzrd") is None, reason="datavzrd is not installed")
def test_render_writes_the_tsv_and_the_site(tmp_path):
    summary = summarize(
        "alfred", {"json": [str(FIXTURES / "lib1.target_coverage.json.gz")]}, ["lib1"]
    )

    qc_summary.render(
        summary, "qc: Target coverage", str(tmp_path / "qc_summary.tsv"), str(tmp_path / "html")
    )

    assert (tmp_path / "html" / "index.html").is_file()
    tsv = pl.read_csv(tmp_path / "qc_summary.tsv", separator="\t")
    assert tsv["library"].to_list() == ["lib1"]
