# -*- coding: utf-8 -*-
"""The fingerprint checks of an ``ngs_data_qc`` task with ``tool: ngs_chew``: the sex and the
pedigree of the sample sheet against the data.

``ngs-chew compare`` writes the relatedness (peddy's estimate), the IBS0 count and the number of sites both libraries cover for every pair of libraries; ``ngs-chew stats`` writes per library the alt-carrying over homozygous ratio at chrX sites (``chrx_het_hom``) and the fractions of reads on chrX and chrY.
The sample sheet facts come as options: ``libraries`` (library, donor, pedigree, sex) and ``donors`` (every donor of the pedigrees, sequenced or not, with its parents), so that the expected relatedness of a pair follows from the pedigree (kinship over the donor graph, unsequenced donors included) and each library's parents can be checked.

Sex is predicted without a fixed cutoff, since the ratio depends on the kit: a library gets the sheet sex whose cohort median ``chrx_het_hom`` is nearer, which needs libraries of both sexes with a known sheet sex; otherwise there is no prediction.
The chrY read fraction is shown next to it as the second evidence.
A pair agrees with the sheet when its relatedness is within :data:`TOLERANCE` of the expected one; a parent and child must also share an allele at nearly every site (IBS0 rate at most :data:`PARENT_MAX_IBS0`), which only errors break; a pair of different pedigrees with a relatedness of at least :data:`UNEXPECTED_RELATEDNESS` is an unexpected relative (a swap, a duplicate, an undeclared relation).
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    MUTED,
    STATUS,
    Chart,
    Column,
    Summary,
    Table,
    bars,
    pills,
    ticks,
)

#: How far the relatedness of a pair may lie from the sheet's expectation
TOLERANCE = 0.15

#: IBS0 rate up to which a parent and its child share an allele at every site but errors
PARENT_MAX_IBS0 = 0.01

#: Relatedness from which a pair of different pedigrees is an unexpected relative
UNEXPECTED_RELATEDNESS = 0.2

#: The outcomes of a check and their colors; the label always carries the meaning
CHECKS = (
    ("ok", STATUS["good"]),
    ("mismatch", STATUS["critical"]),
    ("unrelated", STATUS["critical"]),
    ("no prediction", MUTED),
    ("no sheet sex", MUTED),
    ("not sequenced", MUTED),
    ("none", MUTED),
)

SUMMARY_COLUMNS = {
    "pedigree": Column(),
    "sex": Column(label="sex (sheet)"),
    "predicted sex": Column(),
    "sex check": pills(*zip(*CHECKS)),
    "chrx_het_hom": ticks(label="chrX alt/hom", precision=3),
    "chrx_frac": ticks(label="chrX reads", precision=4),
    "chry_frac": ticks(label="chrY reads", precision=4),
    "father": Column(label="father (sheet)"),
    "father check": pills(*zip(*CHECKS)),
    "mother": Column(label="mother (sheet)"),
    "mother check": pills(*zip(*CHECKS)),
    "relatedness check": pills(*zip(*CHECKS)),
    "unexpected relatives": Column(),
    "hets": Column(display="detail"),
    "hom_alts": Column(display="detail"),
    "hom_refs": Column(display="detail"),
    "var_het": Column(display="detail", precision=4),
}

PAIR_COLUMNS = {
    "pedigree": Column(),
    "a": Column(),
    "b": Column(),
    "expected": ticks(label="relatedness (sheet)", precision=2, domain=(0, 1)),
    "relatedness": ticks(precision=2, domain=(0, 1)),
    "difference": ticks(precision=2),
    "ibs0_rate": ticks(label="IBS0 rate", precision=4),
    "ibs0": bars(label="IBS0 sites", display="available"),
    "sites": bars(display="available"),
    "check": pills(*zip(*CHECKS)),
}


def read_matrix(path: str, value: str) -> pl.DataFrame:
    """The pairs ``a < b`` of a square matrix of ``ngs-chew compare`` with the ``value``."""
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    name = frame.columns[0]
    libraries = frame.columns[1:]
    rows = []
    for row in frame.iter_rows(named=True):
        for other in libraries:
            if row[name] < other:
                rows.append({"a": row[name], "b": other, value: float(row[other])})
    return pl.DataFrame(rows, schema={"a": pl.String, "b": pl.String, value: pl.Float64})


def read_stats(path: str) -> pl.DataFrame:
    """The statistics of ``ngs-chew stats``, the numbers as numbers (``-`` is missing)."""
    frame = pl.read_csv(path, separator="\t", infer_schema=False, null_values=["-"])
    numbers = [c for c in frame.columns if c != "sample_name"]
    return frame.rename({"sample_name": "library"}).with_columns(
        pl.col(c).cast(pl.Float64) for c in numbers
    )


def kinships(donors: list[dict[str, Any]]):
    """A function ``(donor_a, donor_b) -> kinship coefficient`` over the donor graph of one
    pedigree; founders are unrelated, an unknown parent (``"0"``) is a founder of its own."""
    parents = {d["donor"]: (d.get("father") or "0", d.get("mother") or "0") for d in donors}

    @lru_cache(maxsize=None)
    def depth(donor: str) -> int:
        father, mother = parents.get(donor, ("0", "0"))
        return 0 if father == "0" and mother == "0" else 1 + max(depth(father), depth(mother))

    @lru_cache(maxsize=None)
    def phi(a: str, b: str) -> float:
        if a == "0" or b == "0":
            return 0.0
        if a == b:
            father, mother = parents.get(a, ("0", "0"))
            return 0.5 * (1 + phi(father, mother))
        if depth(a) < depth(b):  # expand the deeper one, which cannot be an ancestor
            a, b = b, a
        father, mother = parents.get(a, ("0", "0"))
        return 0.5 * (phi(father, b) + phi(mother, b))

    return phi


def predicted_sexes(stats: pl.DataFrame, sexes: dict[str, str]) -> dict[str, str]:
    """Library -> the sheet sex whose cohort median ``chrx_het_hom`` is nearer; empty without
    libraries of both sexes with a known sheet sex."""
    medians = {}
    for sex in ("female", "male"):
        values = [
            row["chrx_het_hom"]
            for row in stats.iter_rows(named=True)
            if sexes.get(row["library"]) == sex and row["chrx_het_hom"] is not None
        ]
        if values:
            medians[sex] = pl.Series(values).median()
    if len(medians) < 2:
        return {}
    return {
        row["library"]: min(medians, key=lambda sex: abs(row["chrx_het_hom"] - medians[sex]))
        for row in stats.iter_rows(named=True)
        if row["chrx_het_hom"] is not None
    }


def summarize(
    files: dict[str, list[str]],
    labels: list[str],
    libraries: list[dict[str, Any]] = (),
    donors: list[dict[str, Any]] = (),
) -> Summary:
    stats = read_stats(files["stats"][0])
    pairs = (
        read_matrix(files["relatedness"][0], "relatedness")
        .join(read_matrix(files["ibs0"][0], "ibs0"), on=["a", "b"], how="left")
        .join(read_matrix(files["mask"][0], "sites"), on=["a", "b"], how="left")
        .with_columns(
            pl.when(pl.col("sites") > 0)
            .then(pl.col("ibs0") / pl.col("sites"))
            .otherwise(None)
            .alias("ibs0_rate")
        )
    )
    by_library = {entry["library"]: entry for entry in libraries}
    by_donor = {(entry["pedigree"], entry["donor"]): entry for entry in donors}
    sexes = {lib: entry.get("sex", "unknown") for lib, entry in by_library.items()}
    predicted = predicted_sexes(stats, sexes)
    phi_of = {
        pedigree: kinships([d for d in donors if d["pedigree"] == pedigree])
        for pedigree in {d["pedigree"] for d in donors}
    }

    def expected_relatedness(a: str, b: str) -> float | None:
        la, lb = by_library.get(a), by_library.get(b)
        if la is None or lb is None:
            return None
        if la["pedigree"] != lb["pedigree"]:
            return 0.0
        phi = phi_of.get(la["pedigree"])
        return 2 * phi(la["donor"], lb["donor"]) if phi else None

    def pair_check(expected, relatedness, ibs0_rate, parent_child) -> str:
        if expected is None or relatedness is None:
            return "none"
        if abs(relatedness - expected) > TOLERANCE:
            return "mismatch"
        if parent_child and ibs0_rate is not None and ibs0_rate > PARENT_MAX_IBS0:
            return "mismatch"
        return "ok"

    def parent_of(library: str, which: str) -> tuple[str, str]:
        """The parent donor named in the sheet and its library (``""`` if unsequenced)."""
        entry = by_library.get(library)
        donor = by_donor.get((entry["pedigree"], entry["donor"])) if entry else None
        parent = donor.get(which) if donor else "0"
        if not parent or parent == "0":
            return "", ""
        parent_donor = by_donor.get((entry["pedigree"], parent))
        return parent, (parent_donor or {}).get("library") or ""

    pair_rows = []
    for row in pairs.iter_rows(named=True):
        a, b = row["a"], row["b"]
        expected = expected_relatedness(a, b)
        parent_child = any(
            parent_of(child, which)[1] == parent
            for child, parent in ((a, b), (b, a))
            for which in ("father", "mother")
        )
        same = by_library.get(a, {}).get("pedigree") == by_library.get(b, {}).get("pedigree")
        pedigree = by_library.get(a, {}).get("pedigree", "") if same else "across pedigrees"
        pair_rows.append(
            {
                "pedigree": pedigree,
                "a": a,
                "b": b,
                "expected": expected,
                "relatedness": row["relatedness"],
                "difference": None if expected is None else row["relatedness"] - expected,
                "ibs0_rate": row["ibs0_rate"],
                "ibs0": row["ibs0"],
                "sites": row["sites"],
                "check": pair_check(expected, row["relatedness"], row["ibs0_rate"], parent_child),
                "same_pedigree": same,
            }
        )
    pair_frame = pl.DataFrame(
        pair_rows,
        schema={
            "pedigree": pl.String,
            "a": pl.String,
            "b": pl.String,
            "expected": pl.Float64,
            "relatedness": pl.Float64,
            "difference": pl.Float64,
            "ibs0_rate": pl.Float64,
            "ibs0": pl.Float64,
            "sites": pl.Float64,
            "check": pl.String,
            "same_pedigree": pl.Boolean,
        },
    )
    pairs_of = {}
    for row in pair_rows:
        pairs_of.setdefault(row["a"], []).append((row["b"], row))
        pairs_of.setdefault(row["b"], []).append((row["a"], row))

    rows = []
    for row in stats.iter_rows(named=True):
        library = row["library"]
        entry = by_library.get(library, {})
        sheet_sex = entry.get("sex", "unknown")
        prediction = predicted.get(library, "")
        if sheet_sex not in ("female", "male"):
            sex_check = "no sheet sex"
        elif not prediction:
            sex_check = "no prediction"
        else:
            sex_check = "ok" if prediction == sheet_sex else "mismatch"
        checks: dict[str, Any] = {}
        for which in ("father", "mother"):
            parent, parent_library = parent_of(library, which)
            checks[which] = parent
            if not parent:
                checks[f"{which} check"] = "none"
            elif not parent_library:
                checks[f"{which} check"] = "not sequenced"
            else:
                pair = next(
                    (p for other, p in pairs_of.get(library, []) if other == parent_library), None
                )
                checks[f"{which} check"] = (
                    "none" if pair is None else "ok" if pair["check"] == "ok" else "unrelated"
                )
        within = [p for _, p in pairs_of.get(library, []) if p["same_pedigree"]]
        if not within or all(p["check"] == "none" for p in within):
            relatedness_check = "none"
        else:
            relatedness_check = (
                "ok" if all(p["check"] in ("ok", "none") for p in within) else "mismatch"
            )
        unexpected = sorted(
            f"{other} ({p['relatedness']:.2f})"
            for other, p in pairs_of.get(library, [])
            if not p["same_pedigree"]
            and p["relatedness"] is not None
            and p["relatedness"] >= UNEXPECTED_RELATEDNESS
        )
        rows.append(
            {
                "library": library,
                "pedigree": entry.get("pedigree", ""),
                "sex": sheet_sex,
                "predicted sex": prediction,
                "sex check": sex_check,
                "chrx_het_hom": row["chrx_het_hom"],
                "chrx_frac": row["chrx_frac"],
                "chry_frac": row["chry_frac"],
                **checks,
                "relatedness check": relatedness_check,
                "unexpected relatives": ", ".join(unexpected),
                "hets": row["hets"],
                "hom_alts": row["hom_alts"],
                "hom_refs": row["hom_refs"],
                "var_het": row["var_het"],
            }
        )
    summary = pl.DataFrame(
        rows, schema={"library": pl.String, **{c: _dtype(c) for c in SUMMARY_COLUMNS}}
    )

    prediction_note = (
        "Predicted sex is the sheet sex whose cohort median chrX alt/hom ratio is nearer. "
        if predicted
        else "No sex is predicted: the cohort lacks libraries of both sexes with a known sheet "
        "sex, which the prediction compares with. "
    )
    scatter_frame = summary.select("library", "sex", "predicted sex", "chrx_het_hom", "chry_frac")
    shown_pairs = pair_frame.filter(
        pl.col("same_pedigree") | (pl.col("relatedness") >= UNEXPECTED_RELATEDNESS)
    ).drop("same_pedigree")
    return Summary(
        table=Table(
            "libraries",
            "Fingerprint checks",
            summary,
            desc="One row per library: the sheet's sex against the sex predicted from the chrX "
            "alt/hom ratio of the fingerprint (with the chrX and chrY read fractions as further "
            "evidence), the sheet's father and mother against the relatedness and IBS0 rate of "
            "the pairs, whether every pair within the pedigree has the relatedness the sheet "
            f"implies (within {TOLERANCE}), and libraries of other pedigrees with a relatedness "
            f"of at least {UNEXPECTED_RELATEDNESS}. " + prediction_note + "The details hold "
            "the autosomal counts and var(het).",
            columns=SUMMARY_COLUMNS,
        ),
        views=[
            Chart(
                "sex",
                "Sex: chrX alt/hom ratio against chrY reads, colored by the sheet's sex",
                scatter_frame,
                x="chrx_het_hom",
                y="chry_frac",
                mark="scatter",
                series="sex",
                tooltip=("library", "predicted sex"),
                x_title="chrX alt/hom ratio",
                y_title="fraction of reads on chrY",
                link=None,
                desc="One dot per library, colored by the sex of the sample sheet. Males sit "
                "low on the chrX ratio and high on chrY reads, females the other way; a dot of "
                "the wrong color in a cluster is a sex mismatch, a dot between the clusters a "
                "library to look at.",
            ),
            Chart(
                "relatedness",
                "Relatedness of the pairs: sample sheet against fingerprints",
                shown_pairs.with_columns((pl.col("a") + " / " + pl.col("b")).alias("pair")),
                x="expected",
                y="relatedness",
                mark="scatter",
                series="pedigree",
                tooltip=("pair", "check"),
                x_title="relatedness implied by the sample sheet",
                y_title="relatedness from the fingerprints",
                link=None,
                desc="One dot per pair of a pedigree, and the pairs across pedigrees that are "
                "related at all. Pairs on the diagonal are as the sheet says; a dot far off it "
                "is a swapped, unrelated or duplicated library.",
            ),
            Table(
                "pairs",
                "Pairs",
                shown_pairs,
                desc="Every pair of libraries within a pedigree, and the pairs across pedigrees "
                f"with a relatedness of at least {UNEXPECTED_RELATEDNESS}: the relatedness the "
                "sheet implies (kinship over the pedigree, unsequenced donors included), the "
                "relatedness and IBS0 rate of the fingerprints, and the check.",
                columns=PAIR_COLUMNS,
                pinned=("pedigree", "a", "b"),
            ),
        ],
    )


def _dtype(column: str):
    return (
        pl.Float64
        if column
        in {"chrx_het_hom", "chrx_frac", "chry_frac", "hets", "hom_alts", "hom_refs", "var_het"}
        else pl.String
    )
