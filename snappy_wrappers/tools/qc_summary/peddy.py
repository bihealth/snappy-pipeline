# -*- coding: utf-8 -*-
"""The peddy checks of a ``variant_checking`` task, one set of files per pedigree.

The summary table has one row per sample with its pedigree, the sex check (sex of the PED file,
predicted sex, error), the heterozygosity check (depth, heterozygosity ratio, outliers, the
inter-decile range of the B allele frequency, call rate) and the ancestry prediction; the pairs
table lists every pair of a pedigree with the relatedness from the PED file and from the
genotypes; a scatter chart puts the pairs on expected versus observed relatedness, so a pair off
the diagonal (a sample swap, an unknown relation, a duplicate) stands out.
"""

from __future__ import annotations

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    Chart,
    Column,
    Summary,
    Table,
    bars,
    metric_strips,
    numeric,
    ticks,
)

#: Columns of the summary table, in order, with their rendering
SUMMARY_COLUMNS = {
    "pedigree": Column(),
    "ped_sex": Column(label="sex (PED)"),
    "predicted_sex": Column(),
    "sex_error": Column(label="sex error"),
    "het_ratio": ticks(label="heterozygosity ratio"),
    "ratio_outlier": Column(label="heterozygosity outlier"),
    "mean_depth": ticks(label="mean depth", precision=1),
    "depth_outlier": Column(),
    "idr_baf": ticks(label="IDR of the B allele frequency"),
    "call_rate": ticks(),
    "ancestry-prediction": Column(label="ancestry"),
    "ancestry-prob": ticks(label="ancestry probability"),
    "sampled_sites": Column(detail=True),
    "median_depth": Column(detail=True),
    "het_count": Column(detail=True),
    "hom_ref_count": Column(detail=True),
    "hom_alt_count": Column(detail=True),
    "sex_het_ratio": Column(detail=True, label="chrX heterozygosity ratio"),
    "p10": Column(detail=True),
    "p90": Column(detail=True),
    "PC1": Column(detail=True),
    "PC2": Column(detail=True),
    "PC3": Column(detail=True),
    "PC4": Column(detail=True),
}

#: The metrics of the strip chart, with their labels
METRICS = {
    "het_ratio": "heterozygosity ratio",
    "mean_depth": "mean depth",
    "idr_baf": "IDR of the B allele frequency",
    "call_rate": "call rate",
    "sex_het_ratio": "chrX heterozygosity ratio",
}

#: Columns of the pairs table, in order
PAIR_COLUMNS = {
    "pedigree": Column(),
    "sample_a": Column(),
    "sample_b": Column(),
    "pedigree_relatedness": ticks(label="relatedness (PED)", precision=2),
    "rel": ticks(label="relatedness (genotypes)", precision=2),
    "rel_difference": ticks(label="difference", precision=2),
    "pedigree_parents": Column(label="parent and child (PED)"),
    "predicted_parents": Column(label="parent and child (genotypes)"),
    "parent_error": Column(),
    "sample_duplication_error": Column(label="duplication error"),
    "ibs0": bars(),
    "ibs2": bars(),
    "shared_hets": Column(detail=True),
    "hets_a": Column(detail=True),
    "hets_b": Column(detail=True),
    "n": Column(detail=True, label="sites"),
}


def read_csv(path: str) -> pl.DataFrame:
    """A peddy CSV with its numbers as numbers and its booleans as text."""
    return numeric(pl.read_csv(path, separator=",", infer_schema=False))


def summarize(files: dict[str, list[str]], labels: list[str]) -> Summary:
    samples: list[pl.DataFrame] = []
    pairs: list[pl.DataFrame] = []
    for pedigree, ped_check, sex_check, het_check in zip(
        labels, files["ped_check"], files["sex_check"], files["het_check"]
    ):
        sex = read_csv(sex_check).rename({"error": "sex_error", "het_ratio": "sex_het_ratio"})
        het = read_csv(het_check)
        frame = het.join(sex, on="sample_id", how="full", coalesce=True)
        samples.append(
            frame.rename({"sample_id": "library"}).with_columns(pl.lit(pedigree).alias("pedigree"))
        )
        pairs.append(read_csv(ped_check).with_columns(pl.lit(pedigree).alias("pedigree")))

    summary = pl.concat(samples, how="diagonal") if samples else pl.DataFrame({"library": []})
    present = [c for c in SUMMARY_COLUMNS if c in summary.columns]
    summary = summary.select("library", *present)
    pair_frame = (
        pl.concat(pairs, how="diagonal")
        if pairs
        else pl.DataFrame({"pedigree": [], "sample_a": [], "sample_b": []})
    )
    pair_present = [c for c in PAIR_COLUMNS if c in pair_frame.columns]
    pair_frame = pair_frame.select(*pair_present).with_columns(
        (pl.col("sample_a") + " / " + pl.col("sample_b")).alias("pair")
    )
    return Summary(
        table=Table(
            "summary",
            "Pedigree checks",
            summary,
            desc="peddy's checks of each sample: the sex of the PED file against the sex "
            "predicted from chrX, the heterozygosity ratio and depth with their outlier flags, "
            "the inter-decile range of the B allele frequency, the call rate and the ancestry "
            "predicted from the thousand genomes principal components. The details hold the "
            "counts.",
            columns={c: SUMMARY_COLUMNS[c] for c in present},
        ),
        views=[
            metric_strips(
                summary,
                METRICS,
                "libraries",
                "Samples by metric",
                desc="One tick per sample on each metric's range: a tick apart from the others "
                "is a sample to look at; the tooltip names it.",
            ),
            Chart(
                "relatedness",
                "Relatedness of the pairs: PED file versus genotypes",
                pair_frame,
                x="pedigree_relatedness",
                y="rel",
                mark="scatter",
                series="pedigree",
                tooltip=("pair",),
                x_title="relatedness in the PED file",
                y_title="relatedness from the genotypes",
                link=None,
                desc="One dot per pair of samples of a pedigree. Pairs on the diagonal are as "
                "the PED file says; a dot far off it is a swapped, unrelated or duplicated "
                "sample.",
            ),
            Table(
                "pairs",
                "Pairs",
                pair_frame.drop("pair"),
                desc="peddy's `ped_check`: every pair of samples of a pedigree with the "
                "relatedness from the PED file and from the genotypes, whether the pair is "
                "parent and child in both, and the sites with no shared allele (IBS0) and with "
                "the same genotype (IBS2).",
                columns={c: PAIR_COLUMNS[c] for c in pair_present},
                pinned=("pedigree", "sample_a", "sample_b"),
            ),
        ],
    )
