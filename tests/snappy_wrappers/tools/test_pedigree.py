# -*- coding: utf-8 -*-
"""Tests of the pedigree plot of the calls report: the members with their unknown parents,
the generations and positions of the layout, the edges and the plot."""

import json
from pathlib import Path

import fastjsonschema

from snappy_wrappers.tools.calls_report import Options, summarize
from snappy_wrappers.tools.calls_report.pedigree import (
    Member,
    generations,
    layout,
    members_of,
    nodes_and_edges,
    plot,
)
from snappy_wrappers.tools.qc_summary import config

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"

#: The simulated pedigree: grandparents (one unsequenced), parents (the mother's father
#: unknown), the index and a sibling
PEDIGREE = [
    {"name": "pgf-N1", "label": "pgf", "sex": "male", "sample": "pgf-N1"},
    {"name": "pgm", "label": "pgm", "sex": "female"},
    {
        "name": "fa-N1",
        "label": "fa",
        "father": "pgf-N1",
        "mother": "pgm",
        "sex": "male",
        "affected": False,
        "sample": "fa-N1",
    },
    {"name": "mgm-N1", "label": "mgm", "sex": "female", "sample": "mgm-N1"},
    {"name": "mo-N1", "label": "mo", "mother": "mgm-N1", "sex": "female", "sample": "mo-N1"},
    {
        "name": "idx-N1",
        "label": "idx",
        "father": "fa-N1",
        "mother": "mo-N1",
        "sex": "female",
        "affected": True,
        "sample": "idx-N1",
    },
    {
        "name": "sib1-N1",
        "label": "sib1",
        "father": "fa-N1",
        "mother": "mo-N1",
        "sex": "male",
        "affected": False,
        "sample": "sib1-N1",
    },
]
SAMPLES = ["pgf-N1", "fa-N1", "mgm-N1", "mo-N1", "idx-N1", "sib1-N1"]


def test_members_get_an_unknown_parent_and_lose_samples_the_vcf_lacks():
    members = members_of(PEDIGREE, [s for s in SAMPLES if s != "sib1-N1"])

    by_name = {member.name: member for member in members}
    assert by_name["mo-N1"].father == "mo-N1.father"
    assert by_name["mo-N1.father"] == Member(name="mo-N1.father", label="?", sex="male")
    assert by_name["pgm"].sample is None and by_name["fa-N1"].sample == "fa-N1"
    assert by_name["sib1-N1"].sample is None  # not in the VCF
    assert by_name["idx-N1"].affected is True and by_name["pgf-N1"].affected is None


def test_generations_put_partners_in_one_row_and_children_below():
    members = members_of(PEDIGREE, SAMPLES)

    rows = generations(members)

    assert rows["pgf-N1"] == rows["pgm"] == rows["mgm-N1"] == rows["mo-N1.father"] == 0
    assert rows["fa-N1"] == rows["mo-N1"] == 1
    assert rows["idx-N1"] == rows["sib1-N1"] == 2


def test_layout_centers_children_under_their_parents_and_keeps_nodes_apart():
    members = members_of(PEDIGREE, SAMPLES)

    positions = layout(members)

    x = {name: value[0] for name, value in positions.items()}
    assert x["pgf-N1"] + 1 == x["pgm"]  # the father left of the mother
    assert x["fa-N1"] == (x["pgf-N1"] + x["pgm"]) / 2  # the child under its parents
    assert (
        x["mo-N1"] == (x["mgm-N1"] + x["mo-N1.father"]) / 2
    )  # under her mother and unknown father
    assert {x["idx-N1"], x["sib1-N1"]} == {
        (x["fa-N1"] + x["mo-N1"]) / 2 - 0.5,
        (x["fa-N1"] + x["mo-N1"]) / 2 + 0.5,
    }
    assert min(x.values()) == 0
    for row in (0, 1, 2):
        values = sorted(v for name, (v, g) in positions.items() if g == row)
        assert all(b - a >= 1 for a, b in zip(values, values[1:]))
    # A founder partner sits next to its partner, the father to the left
    trio = layout(
        members_of(
            PEDIGREE[2:3]
            + PEDIGREE[5:6]
            + [{"name": "mo-N1", "label": "mo", "sex": "female", "sample": "mo-N1"}],
            ["fa-N1", "mo-N1", "idx-N1"],
        )
    )
    assert trio["fa-N1"][0] + 1 == trio["mo-N1"][0] and trio["idx-N1"][0] == trio["fa-N1"][0] + 0.5


def test_nodes_and_edges_carry_the_glyphs_and_link_couples_to_children():
    members = members_of(PEDIGREE, SAMPLES)

    nodes, edges = nodes_and_edges(members)

    by_name = {node["name"]: node for node in nodes}
    assert by_name["idx-N1"]["sex"] == "female" and by_name["idx-N1"]["status"] == "affected"
    assert by_name["pgm"]["observed"] is False and by_name["pgm"]["status"] == "unknown"
    # The stroke carries the sex, the fill the status: an affected member in its sex color
    assert by_name["idx-N1"]["stroke_color"] == by_name["idx-N1"]["fill_color"] == "#eb6834"
    assert (
        by_name["fa-N1"]["stroke_color"] == "#1baf7a"
        and by_name["fa-N1"]["fill_color"] == "#ffffff"
    )
    assert by_name["pgm"]["fill_color"] == "#f0efec"
    assert by_name["fa-N1"]["y"] == -1 and by_name["pgf-N1"]["y"] == 0
    # Orthogonal segments: one per couple between the partners, a drop to the sibship line,
    # the sibship line of the two children, and a drop to each child
    from snappy_wrappers.tools.calls_report.pedigree import SIBSHIP

    partner_lines = [e for e in edges if e["y"] == e["y2"] and e["y"] % 1 == 0]
    child_drops = [e for e in edges if e["x"] == e["x2"] and e["y"] % 1 != 0]
    sibship_lines = [e for e in edges if e["y"] == e["y2"] and e["y"] % 1 != 0]
    assert len(partner_lines) == 3  # pgf-pgm, mo.father-mgm, fa-mo
    assert len(child_drops) == 4  # fa, mo, idx, sib1
    assert len(sibship_lines) == 1  # idx and sib1 share one; the single children need none
    assert all(e["y"] == e["y2"] + 1 - SIBSHIP for e in child_drops)


def test_the_plot_is_a_layered_spec_whose_data_function_reads_the_sample_columns():
    members = members_of(PEDIGREE, SAMPLES)

    custom = plot(members)

    spec, data = custom["spec"], custom["data"]
    assert spec["layer"][0]["mark"]["type"] == "rule"  # edges first, below the nodes
    assert spec["width"] > 0 and spec["height"] > 0
    assert '"sample": "idx-N1"' in data and 'row[node.sample + ": AF"]' in data
    assert 'row[node.sample + ": observations"]' in data  # the reads by strength
    assert data.count('"kind": "legend"') == 6 and '"kind": "note"' in data  # the legend
    assert spec["width"] >= 3.3 * 96  # room for the legend's three entries per row
    assert [layer["mark"]["type"] for layer in spec["layer"]].count(
        "rect"
    ) == 3  # AF bars and reads
    assert "data-value" in data  # values that datavzrd renders as HTML still parse
    assert 'row["P(" + event + ")"]' in data


def test_the_calls_report_adds_the_pedigree_column_in_the_row_details():
    pedigree = [
        {
            "name": "case001subregion-N1-DNA1-WES1",
            "label": "case001subregion",
            "father": "case001subregionFather-N1-DNA1-WES1",
            "mother": "case001subregionMother-N1-DNA1-WES1",
            "sex": "male",
            "affected": True,
            "sample": "case001subregion-N1-DNA1-WES1",
        },
        {
            "name": "case001subregionFather-N1-DNA1-WES1",
            "label": "case001subregionFather",
            "sex": "male",
            "affected": False,
            "sample": "case001subregionFather-N1-DNA1-WES1",
        },
        {
            "name": "case001subregionMother-N1-DNA1-WES1",
            "label": "case001subregionMother",
            "sex": "female",
            "affected": False,
            "sample": "case001subregionMother-N1-DNA1-WES1",
        },
    ]

    summary = summarize(str(FIXTURES / "calls.vcf"), Options(pedigree=tuple(pedigree)))
    table = summary.table

    assert table.columns["pedigree"].display == "detail"
    assert "pedigree" in table.exclude  # a custom plot: in the event views only
    order = list(table.columns)
    assert order.index("pedigree") == order.index("P(absent)") + 1  # first of the row details
    assert order[-1] == "hints"  # the hints last
    # The per-sample columns step back behind the pedigree
    assert table.columns["case001subregion-N1-DNA1-WES1: AF"].display == "available"
    assert table.columns["case001subregion-N1-DNA1-WES1: GT"].display == "available"
    assert (
        summarize(str(FIXTURES / "calls.vcf"))
        .table.columns["case001subregion-N1-DNA1-WES1: AF"]
        .display
        is None
    )
    assert "`pedigree` (row details, first)" in table.desc
    present = next(view for view in summary.views if view.name == "present")
    cfg = config([table, *summary.views], "calls", "/data", summary.default_view)
    assert "pedigree" in cfg["views"]["present"]["render-table"]["columns"]
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)
    # Without a pedigree, or with one member, there is no column
    assert "pedigree" not in summarize(str(FIXTURES / "calls.vcf")).table.frame.columns
    assert (
        "pedigree"
        not in summarize(
            str(FIXTURES / "calls.vcf"), Options(pedigree=(pedigree[0],))
        ).table.frame.columns
    )
    assert present.frame["pedigree"].to_list() == [""] * 6
