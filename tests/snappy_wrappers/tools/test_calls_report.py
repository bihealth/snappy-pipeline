# -*- coding: utf-8 -*-
"""Tests of the calls report: the VCF reader, the varlociraptor and the generic sample columns,
the compact cells and the datavzrd configuration."""

import json
from pathlib import Path

import fastjsonschema
import pytest

from snappy_wrappers.tools.calls_report import Options, summarize, variant_type
from snappy_wrappers.tools.calls_report.varlociraptor import (
    biases,
    compact_afd,
    observations,
    phred_to_probability,
)
from snappy_wrappers.tools.calls_report.vcf import header_attributes, read_vcf
from snappy_wrappers.tools.qc_summary import AGAINST, CELL_CONFIG, SEQUENTIAL, VEGA_LITE, config

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
CALLS = FIXTURES / "calls.vcf"

#: A VCF of another caller: a multi-allelic record and a record without data for one sample
GENERIC_VCF = """##fileformat=VCFv4.2
##contig=<ID=1>
##FILTER=<ID=low_qual,Description="low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths, comma-separated">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\trs1\tG\tT,C\t50.5\tPASS\tDP=40\tGT:AD:DP:GQ\t1/2:10,20,10:40:99\t0/1:30,10,0:40:60
1\t200\t.\tA\tAT\t.\tlow_qual\tDP=10\tGT:AD:DP:GQ\t0/0:10,0:10:30\t./.:.:.:.
"""


def _validate(summary):
    """The datavzrd configuration of ``summary`` is valid against the schema."""
    cfg = config([summary.table, *summary.views], "calls", "/data", summary.default_view)
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)
    return cfg


def test_read_vcf_takes_the_header_fields_and_one_record_per_alt(tmp_path):
    vcf = read_vcf(str(CALLS))

    assert vcf.samples == [
        "case001subregion-N1-DNA1-WES1",
        "case001subregionFather-N1-DNA1-WES1",
        "case001subregionMother-N1-DNA1-WES1",
    ]
    assert [name for name in vcf.info if name.startswith("PROB_")] == [
        "PROB_DENOVO",
        "PROB_PRESENT",
        "PROB_ARTIFACT",
        "PROB_ABSENT",
    ]
    assert "Description" in vcf.info["PROB_DENOVO"] and "AF" in vcf.formats
    assert len(vcf.records) == 6
    record = vcf.records[0]
    assert (record.chromosome, record.position, record.alt_index, record.alts) == (
        "12",
        23687354,
        0,
        1,
    )
    assert record.info["PROB_PRESENT"] == "0"
    assert record.samples[vcf.samples[0]]["GT"] == "0/1"

    path = tmp_path / "generic.vcf"
    path.write_text(GENERIC_VCF, encoding="utf-8")
    generic = read_vcf(str(path))

    assert generic.contigs == ["1"] and list(generic.formats) == ["GT", "AD", "DP", "GQ"]
    assert [(r.alt, r.alt_index, r.alts) for r in generic.records] == [
        ("T", 0, 2),
        ("C", 1, 2),
        ("AT", 0, 1),
    ]
    assert header_attributes('ID=AD,Number=R,Description="Allelic depths, comma-separated"') == {
        "ID": "AD",
        "Number": "R",
        "Description": "Allelic depths, comma-separated",
    }


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


def test_varlociraptor_calls_get_the_probabilities_and_the_per_sample_columns():
    summary = summarize(str(CALLS), Options(page_size=10, max_in_memory_rows=1000))
    table = summary.table
    frame = table.frame

    assert table.name == "all" and table.pinned == ("chromosome", "position")
    assert (table.page_size, table.max_in_memory_rows) == (10, 1000)
    assert frame.columns[:9] == [
        "chromosome",
        "position",
        "type",
        "ref",
        "alt",
        "P(denovo)",
        "P(present)",
        "P(artifact)",
        "P(absent)",
    ]
    assert frame.height == 6
    row = frame.row(0, named=True)
    assert (row["chromosome"], row["position"], row["type"]) == ("12", 23687354, "SNV")
    assert row["P(present)"] == 1.0  # PHRED 0
    assert row["P(artifact)"] == 0.0  # PHRED inf
    assert row["P(denovo)"] == pytest.approx(10 ** (-2796.69 / 10))
    index = "case001subregion-N1-DNA1-WES1"
    assert (row[f"{index}: AF"], row[f"{index}: DP"], row[f"{index}: GT"]) == (0.5, 57, "0/1")
    assert row[f"{index}: observations"] == "26V2B|27V2B"
    assert row[f"{index}: AFD"] == "0=924.7,0.5=0.0,1=915.3"
    # The VCF fills neither ID, QUAL nor FILTER, so there are no such columns
    assert not {"id", "quality", "filter"} & set(frame.columns)
    assert table.columns["hints"].display == "detail"

    assert table.columns["P(present)"].plot["heatmap"]["range"] == list(SEQUENTIAL)
    assert table.columns["P(artifact)"].plot["heatmap"]["range"] == list(AGAINST)
    sample = "case001subregionFather-N1-DNA1-WES1"
    assert table.columns[f"{sample}: AF"].plot["ticks"]["domain"] == [0, 1]
    assert table.columns[f"{sample}: DP"].plot["bars"]
    plot = table.columns[f"{sample}: observations"].custom_plot
    assert plot["spec"]["mark"]["type"] == "bar"
    assert plot["spec"]["$schema"] == VEGA_LITE and plot["spec"]["config"] == CELL_CONFIG
    assert 'value.split("|")' in plot["data"]
    # Every strength has its row and the cells of the column share the x scale
    encoding = plot["spec"]["encoding"]
    assert encoding["y"]["scale"]["domain"] == [
        "very strong",
        "strong",
        "positive",
        "barely",
        "equal",
    ]
    most = max(
        record["count"]
        for value in frame[f"{sample}: observations"].to_list()
        for record in observations(*value.split("|"))
    )
    assert encoding["x"]["scale"]["domain"] == [-most, most]
    assert table.columns[f"{sample}: observations"].display == "detail"
    assert table.columns[f"{sample}: AFD"].custom_plot["spec"]["mark"]["type"] == "bar"
    assert table.columns[f"{sample}: AFD"].display == "detail"
    assert table.columns[f"{sample}: biases"].display == "available"

    cfg = _validate(summary)
    column = cfg["views"]["all"]["render-table"]["columns"][f"{sample}: observations"]
    assert json.loads(column["custom-plot"]["spec"])["width"] == 110
    assert cfg["views"]["all"]["max-in-memory-rows"] == 1000


def test_other_callers_get_genotype_frequency_and_depth_columns(tmp_path):
    path = tmp_path / "calls.vcf"
    path.write_text(GENERIC_VCF, encoding="utf-8")

    summary = summarize(str(path))
    frame = summary.table.frame

    assert frame.columns == [
        "chromosome",
        "position",
        "type",
        "ref",
        "alt",
        "S1: GT",
        "S1: AF",
        "S1: alt reads",
        "S1: DP",
        "S1: GQ",
        "S2: GT",
        "S2: AF",
        "S2: alt reads",
        "S2: DP",
        "S2: GQ",
        "id",
        "quality",
        "filter",
    ]
    assert frame["type"].to_list() == ["SNV", "SNV", "INS"]
    assert frame["S1: AF"].to_list() == pytest.approx([0.5, 0.25, 0.0])
    assert frame["S1: alt reads"].to_list() == [20, 10, 0]
    assert frame["S2: AF"].to_list()[:2] == pytest.approx([0.25, 0.0])
    assert frame["S2: AF"][2] is None and frame["S2: DP"][2] is None
    assert frame["id"].to_list() == ["rs1", "rs1", ""]
    assert frame["quality"].to_list() == [50.5, 50.5, None]
    assert frame["filter"].to_list() == ["PASS", "PASS", "low_qual"]
    assert summary.table.columns["filter"].display is None  # a filter besides PASS is shown
    assert summary.table.columns["id"].display == "available"
    assert summary.table.columns["S1: GQ"].display == "available"
    assert summary.table.columns["S1: AF"].plot["ticks"]["domain"] == [0, 1]
    _validate(summary)


def test_summary_of_an_empty_vcf_has_the_columns_and_no_rows(tmp_path):
    header = "\n".join(
        line for line in CALLS.read_text(encoding="utf-8").splitlines() if line.startswith("#")
    )
    path = tmp_path / "empty.vcf"
    path.write_text(header + "\n", encoding="utf-8")

    summary = summarize(str(path))

    assert summary.table.frame.height == 0
    assert "P(present)" in summary.table.frame.columns
    assert summary.table.frame.columns[:5] == ["chromosome", "position", "type", "ref", "alt"]
    _validate(summary)
