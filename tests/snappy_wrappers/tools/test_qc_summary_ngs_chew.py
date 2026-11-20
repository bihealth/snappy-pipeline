# -*- coding: utf-8 -*-
"""Tests of the fingerprint checks: the kinship over the pedigree, the sex prediction, the
checks of the pairs and the parents, and the datavzrd configuration."""

import json
from pathlib import Path

import fastjsonschema
import pytest

from snappy_wrappers.tools.qc_summary import Chart, config, summarize
from snappy_wrappers.tools.qc_summary.ngs_chew import kinships, predicted_sexes, read_stats

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
FILES = {
    kind: [str(FIXTURES / f"cohort.{kind}.txt")]
    for kind in ("relatedness", "ibs0", "mask", "stats")
}


def _lib(short):
    return f"{short}-N1-DNA1-WES1"


def _entry(short, pedigree, sex):
    return {"library": _lib(short), "donor": short, "pedigree": pedigree, "sex": sex}


def _donor(donor, pedigree, father="0", mother="0", sex="unknown", sequenced=True):
    return {
        "donor": donor,
        "pedigree": pedigree,
        "father": father,
        "mother": mother,
        "sex": sex,
        "library": _lib(donor) if sequenced else "",
    }


#: Two pedigrees, a singleton of unknown sex and a duplicate of pB1 in a pedigree of its own;
#: FAM_A has a grandfather through an unsequenced grandmother, FAM_B two declared siblings
LIBRARIES = [
    _entry("gfA", "FAM_A", "male"),
    _entry("faA", "FAM_A", "male"),
    _entry("moA", "FAM_A", "female"),
    _entry("idxA", "FAM_A", "male"),
    _entry("pB1", "FAM_B", "female"),
    _entry("pB2", "FAM_B", "male"),  # female by the data: a sex swap
    _entry("cC1", "FAM_C", "unknown"),
    _entry("dup", "FAM_D", "female"),
]
DONORS = [
    _donor("gfA", "FAM_A", sex="male"),
    _donor("gmA", "FAM_A", sex="female", sequenced=False),
    _donor("faA", "FAM_A", "gfA", "gmA", "male"),
    _donor("moA", "FAM_A", sex="female"),
    _donor("idxA", "FAM_A", "faA", "moA", "male"),
    _donor("fB", "FAM_B", sex="male", sequenced=False),
    _donor("mB", "FAM_B", sex="female", sequenced=False),
    _donor("pB1", "FAM_B", "fB", "mB", "female"),
    _donor("pB2", "FAM_B", "fB", "mB", "male"),
    _donor("cC1", "FAM_C"),
    _donor("dup", "FAM_D", sex="female"),
]


def _summary(libraries=LIBRARIES, donors=DONORS):
    return summarize("ngs_chew", FILES, ["cohort"], libraries=libraries, donors=donors)


def test_kinship_follows_the_pedigree_through_unsequenced_donors():
    phi = kinships([d for d in DONORS if d["pedigree"] == "FAM_A"])

    assert 2 * phi("faA", "idxA") == 0.5  # parent and child
    assert 2 * phi("gfA", "idxA") == 0.25  # grandparent, through the unsequenced grandmother
    assert 2 * phi("gfA", "moA") == 0.0  # unrelated founders
    assert phi("idxA", "idxA") == 0.5  # no inbreeding
    siblings = kinships(
        [_donor("f", "P"), _donor("m", "P"), _donor("a", "P", "f", "m"), _donor("b", "P", "f", "m")]
    )
    assert 2 * siblings("a", "b") == 0.5
    half = kinships([_donor("f", "P"), _donor("a", "P", "f", "m1"), _donor("b", "P", "f", "m2")])
    assert 2 * half("a", "b") == 0.25  # half-siblings, the mothers unknown founders


def test_sex_is_predicted_from_the_nearer_median_of_the_sheets_sexes():
    stats = read_stats(FILES["stats"][0])
    sexes = {entry["library"]: entry["sex"] for entry in LIBRARIES}

    predicted = predicted_sexes(stats, sexes)

    assert predicted[_lib("pB2")] == "female"  # the sheet says male
    assert predicted[_lib("cC1")] == "male"  # no sheet sex, still predicted
    assert predicted[_lib("moA")] == "female" and predicted[_lib("gfA")] == "male"
    assert predicted_sexes(stats, {lib: "male" for lib in sexes}) == {}


def test_summary_checks_sex_parents_and_relatedness_against_the_sheet():
    summary = _summary()
    frame = summary.table.frame
    rows = {row["library"]: row for row in frame.iter_rows(named=True)}

    assert frame.columns[:5] == ["library", "pedigree", "sex", "predicted sex", "sex check"]
    assert rows[_lib("pB2")]["sex check"] == "mismatch"
    assert rows[_lib("cC1")]["sex check"] == "no sheet sex"
    assert rows[_lib("idxA")]["sex check"] == "ok"
    assert (rows[_lib("idxA")]["father"], rows[_lib("idxA")]["father check"]) == ("faA", "ok")
    assert (rows[_lib("idxA")]["mother"], rows[_lib("idxA")]["mother check"]) == ("moA", "ok")
    assert rows[_lib("faA")]["mother check"] == "not sequenced"
    assert rows[_lib("gfA")]["father check"] == "none"
    assert rows[_lib("idxA")]["relatedness check"] == "ok"  # incl. the grandfather at 0.26
    assert rows[_lib("pB1")]["relatedness check"] == "mismatch"  # declared siblings, unrelated
    assert rows[_lib("pB1")]["unexpected relatives"] == f"{_lib('dup')} (0.99)"
    assert rows[_lib("dup")]["unexpected relatives"] == f"{_lib('pB1')} (0.99)"
    assert rows[_lib("cC1")]["relatedness check"] == "none"
    assert rows[_lib("gfA")]["chry_frac"] == pytest.approx(0.0021)
    assert rows[_lib("gfA")]["hets"] == 9000
    assert "Predicted sex is the sheet sex whose cohort median" in summary.table.desc


def test_pairs_table_and_charts():
    summary = _summary()
    views = {view.name: view for view in summary.views}

    pairs = views["pairs"].frame
    assert pairs.height == 8  # the six pairs of FAM_A, the FAM_B pair and the duplicate
    duplicate = pairs.filter(pairs["pedigree"] == "across pedigrees").row(0, named=True)
    assert (duplicate["a"], duplicate["b"], duplicate["check"]) == (
        _lib("dup"),
        _lib("pB1"),
        "mismatch",
    )
    grandparent = pairs.filter((pairs["a"] == _lib("gfA")) & (pairs["b"] == _lib("idxA"))).row(
        0, named=True
    )
    assert (grandparent["expected"], grandparent["relatedness"], grandparent["check"]) == (
        0.25,
        0.26,
        "ok",
    )
    assert grandparent["ibs0_rate"] == pytest.approx(900 / 20000)
    unrelated = pairs.filter((pairs["a"] == _lib("gfA")) & (pairs["b"] == _lib("moA"))).row(
        0, named=True
    )
    assert unrelated["relatedness"] == pytest.approx(-0.0103) and unrelated["check"] == "ok"
    # peddy's estimate leaves [0, 1], so the relatedness ticks have no fixed domain
    assert "domain" not in views["pairs"].columns["relatedness"].plot["ticks"]
    assert views["pairs"].columns["expected"].plot["ticks"]["domain"] == [0, 1]
    assert isinstance(views["sex"], Chart) and views["sex"].series == "sex"
    assert (views["sex"].x, views["sex"].y) == ("chrx_het_hom", "chry_frac")
    assert views["relatedness"].mark == "scatter" and views["relatedness"].x == "expected"
    cfg = config([summary.table, *summary.views], "Fingerprint checks", "/data")
    assert cfg["views"]["libraries"]["render-table"]["columns"]["sex check"]["plot"]["pills"]
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)


def test_a_parent_whose_pair_shares_no_alleles_is_unrelated(tmp_path):
    # Declare pB1 the mother of pB2: the pair has IBS0 at 9 % of the sites and relatedness 0.01
    donors = [d for d in DONORS if d["donor"] != "pB2"] + [
        _donor("pB2", "FAM_B", "fB", "pB1", "male")
    ]

    rows = {
        row["library"]: row for row in _summary(donors=donors).table.frame.iter_rows(named=True)
    }

    assert (rows[_lib("pB2")]["mother"], rows[_lib("pB2")]["mother check"]) == ("pB1", "unrelated")


def test_without_sheet_facts_the_checks_have_no_expectation():
    summary = _summary(libraries=[], donors=[])
    frame = summary.table.frame

    assert frame["sex check"].unique().to_list() == ["no sheet sex"]
    assert frame["relatedness check"].unique().to_list() == ["none"]
    assert frame["predicted sex"].unique().to_list() == [""]
    assert "No sex is predicted" in summary.table.desc
