# -*- coding: utf-8 -*-
"""Tests of the varlociraptor calls report: the VCF parser, the observation and density records
of the in-cell plots and the datavzrd configuration."""

import json
from pathlib import Path

import fastjsonschema
import polars as pl
import pytest

from snappy_wrappers.tools.qc_summary import AGAINST, CELL_CONFIG, SEQUENTIAL, VEGA_LITE, config
from snappy_wrappers.tools.qc_summary.varlociraptor_calls import (
    biases,
    compact_afd,
    observations,
    phred_to_probability,
    read_calls,
    summarize,
    variant_type,
)

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
CALLS = FIXTURES / "calls.vcf"


def test_read_calls_takes_the_events_from_the_header_and_one_row_per_record():
    calls = read_calls(str(CALLS))

    assert calls.samples == [
        "case001subregion-N1-DNA1-WES1",
        "case001subregionFather-N1-DNA1-WES1",
        "case001subregionMother-N1-DNA1-WES1",
    ]
    assert calls.events == ["DENOVO", "PRESENT", "ARTIFACT", "ABSENT"]
    assert len(calls.rows) == 6
    row = calls.rows[0]
    assert (row["chromosome"], row["position"], row["type"]) == ("12", 23687354, "SNV")
    assert row["prob: present"] == 1.0  # PHRED 0
    assert row["prob: artifact"] == 0.0  # PHRED inf
    assert row["prob: denovo"] == pytest.approx(10 ** (-2796.69 / 10))
    index = "case001subregion-N1-DNA1-WES1"
    assert row[f"{index}: allele frequency"] == 0.5
    assert row[f"{index}: depth"] == 57
    assert row[f"{index}: genotype"] == "0/1"
    assert row[f"{index}: observations"] == "26V2B|27V2B"
    assert row[f"{index}: allele frequency density"] == "0=924.7,0.5=0.0,1=915.3"


def test_observations_split_alt_and_ref_reads_by_strength():
    records = observations("26V2B", "27V2b")

    assert records == [
        {
            "allele": "alt",
            "strength": "very strong",
            "mapping": "certain",
            "count": 26,
            "reads": 26,
        },
        {"allele": "alt", "strength": "barely", "mapping": "certain", "count": 2, "reads": 2},
        {
            "allele": "ref",
            "strength": "very strong",
            "mapping": "certain",
            "count": 27,
            "reads": -27,
        },
        {"allele": "ref", "strength": "barely", "mapping": "uncertain", "count": 2, "reads": -2},
    ]
    assert observations("", ".") == []


def test_compact_afd_thins_to_the_points_and_keeps_the_mode():
    assert compact_afd("0.000=924.70,0.500=-0.00,1.000=915.34") == "0=924.7,0.5=0.0,1=915.3"
    assert compact_afd("0.000=inf,0.500=0.00") == "0.5=0.0"
    assert compact_afd(".") == ""

    dense = ",".join(f"{i / 100:.3f}={abs(37 - i) * 3:.2f}" for i in range(101))
    entries = compact_afd(dense, 10).split(",")

    assert len(entries) == 10
    assert entries[0] == "0=111.0" and entries[-1] == "1=189.0"
    assert "0.37=0.0" in entries  # the mode survives the thinning


def test_biases_name_the_biases_found():
    assert biases({"SB": "+", "ROB": ".", "RPB": ".", "SCB": "$", "HE": ".", "ALB": "."}) == (
        "strand +, softclip $"
    )
    assert biases({"SB": ".", "ROB": "."}) == ""


def test_phred_to_probability_and_variant_type():
    assert phred_to_probability("0") == 1.0
    assert phred_to_probability("10") == pytest.approx(0.1)
    assert phred_to_probability("inf") == 0.0
    assert phred_to_probability(".") is None
    assert variant_type("A", "AT", {}) == "INS"
    assert variant_type("AT", "A", {}) == "DEL"
    assert variant_type("AT", "GC", {}) == "MNV"
    assert variant_type("N", "<DEL>", {"SVTYPE": "DEL"}) == "DEL"


def test_summary_table_has_the_probabilities_and_the_per_sample_columns_with_cell_plots():
    summary = summarize(str(CALLS))
    table = summary.table

    assert table.name == "calls" and table.pinned == ("chromosome", "position")
    assert table.frame.columns[:5] == ["chromosome", "position", "type", "ref", "alt"]
    assert table.frame.height == 6
    assert table.columns["prob: present"].plot["heatmap"]["range"] == list(SEQUENTIAL)
    assert table.columns["prob: artifact"].plot["heatmap"]["range"] == list(AGAINST)
    sample = "case001subregionFather-N1-DNA1-WES1"
    assert table.columns[f"{sample}: allele frequency"].plot["ticks"]["domain"] == [0, 1]
    assert table.columns[f"{sample}: depth"].plot["bars"]
    observation_plot = table.columns[f"{sample}: observations"].custom_plot
    assert observation_plot["spec"]["mark"]["type"] == "bar"
    assert observation_plot["spec"]["$schema"] == VEGA_LITE
    assert observation_plot["spec"]["config"] == CELL_CONFIG
    assert 'value.split("|")' in observation_plot["data"]
    # Every strength has its row and the cells of the column share the x scale
    encoding = observation_plot["spec"]["encoding"]
    assert encoding["y"]["scale"]["domain"] == [
        "very strong",
        "strong",
        "positive",
        "barely",
        "equal",
    ]
    most = max(
        record["count"]
        for value in table.frame[f"{sample}: observations"].to_list()
        for record in observations(*value.split("|"))
    )
    assert encoding["x"]["scale"]["domain"] == [-most, most]
    assert (
        table.columns[f"{sample}: allele frequency density"].custom_plot["spec"]["mark"]["type"]
        == "bar"
    )
    assert f"{sample}: raw observations" not in table.frame.columns
    assert table.frame[f"{sample}: observations"][0] == "26V2B|27V2B"

    cfg = config([table], "calls", "/data")
    column = cfg["views"]["calls"]["render-table"]["columns"][f"{sample}: observations"]
    assert json.loads(column["custom-plot"]["spec"])["width"] == 110
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)


def test_summary_of_an_empty_vcf_has_the_columns_and_no_rows(tmp_path):
    header = "\n".join(
        line for line in CALLS.read_text(encoding="utf-8").splitlines() if line.startswith("#")
    )
    path = tmp_path / "empty.vcf"
    path.write_text(header + "\n", encoding="utf-8")

    summary = summarize(str(path))

    assert summary.table.frame.height == 0
    assert "prob: present" in summary.table.frame.columns
    assert isinstance(summary.table.frame, pl.DataFrame)
