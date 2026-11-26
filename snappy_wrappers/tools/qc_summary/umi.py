# -*- coding: utf-8 -*-
"""The QC files of an ``umi_mapping`` task: the UMI metrics of fgumi and ``samtools`` of the consensus BAM.

The summary table has one row per library: the ``samtools stats`` columns of the consensus BAM (:mod:`.samtools`), then the UMI grouping, the consensus calling, the filter and the yield of the full data.
The charts show the family sizes of the consensus calling, of the grouping and of the position groups, for duplex consensus the AB by BA family sizes of each library, and the yield over the downsampled fractions.
Tables list the UMI counts and, where the task ran ``fgumi correct`` or ``fgumi clip``, their metrics; the samtools views follow.
The option ``consensus`` (``duplex`` or ``simplex``) names the consensus caller, which sets the kinds of families and the yield columns.

The kinds of files and their columns, from the fgumi 0.9.0 sources:

- ``grouping_metrics`` (``fgumi group --metrics``, ``<prefix>.grouping_metrics.txt``): ``accepted_sam_records``, ``discarded_non_pf``, ``discarded_poor_alignment``, ``discarded_ns_in_umi``, ``discarded_umis_to_short`` (``crates/fgumi-metrics/src/group.rs``, ``UmiGroupingMetrics``).
- ``group_family_sizes`` (``<prefix>.family_sizes.txt``): ``family_size``, ``count``, ``fraction``, ``fraction_gt_or_eq_family_size`` (``group.rs``, ``FamilySizeMetrics``).
- ``position_group_sizes`` (``<prefix>.position_group_sizes.txt``): ``position_group_size``, ``count``, ``fraction``, ``fraction_gt_or_eq_position_group_size`` (``group.rs``, ``PositionGroupSizeMetrics``).
- ``family_sizes`` (``fgumi duplex|simplex --metrics``, ``<prefix>.family_sizes.txt``): ``family_size``, then ``count``, ``fraction`` and ``fraction_gt_or_eq_size`` prefixed with ``cs_``, ``ss_`` and, for duplex, ``ds_`` (``crates/fgumi-metrics/src/duplex.rs``, ``FamilySizeMetric``; ``simplex.rs``, ``SimplexFamilySizeMetric``).
- ``duplex_family_sizes`` (``<prefix>.duplex_family_sizes.txt``): ``ab_size``, ``ba_size``, ``count``, ``fraction``, ``fraction_gt_or_eq_size`` (``duplex.rs``, ``DuplexFamilySizeMetric``).
- ``yield_metrics`` (``<prefix>.duplex_yield_metrics.txt``): ``fraction``, ``read_pairs``, ``cs_families``, ``ss_families``, ``ds_families``, ``ds_duplexes``, ``ds_fraction_duplexes``, ``ds_fraction_duplexes_ideal`` (``duplex.rs``, ``DuplexYieldMetric``).
- ``yield_metrics`` (``<prefix>.simplex_yield_metrics.txt``): ``fraction``, ``read_pairs``, ``cs_families``, ``ss_families``, ``mean_ss_family_size``, ``ss_singletons``, ``ss_singleton_fraction``, ``ss_consensus_families`` (``simplex.rs``, ``SimplexYieldMetric``).
- ``umi_counts`` (``<prefix>.umi_counts.txt``): ``umi``, ``raw_observations``, ``raw_observations_with_errors``, ``unique_observations``, ``fraction_raw_observations``, ``fraction_unique_observations`` (``crates/fgumi-metrics/src/shared.rs``, ``UmiMetric``).
- ``consensus_stats`` (``fgumi duplex|simplex --stats``): rows of ``key``, ``value`` and ``description``; the keys are ``raw_reads_considered``, ``raw_reads_rejected``, ``raw_reads_used``, ``frac_raw_reads_used``, one ``raw_reads_rejected_for_<reason>`` per reason and ``consensus_reads_emitted`` (``crates/fgumi-metrics/src/consensus.rs``, ``ConsensusMetrics::to_kv_metrics``).
- ``filter_stats`` (``fgumi filter --stats``): one row of ``total_reads``, ``passed_reads``, ``failed_reads``, ``pass_rate`` (``crates/fgumi-metrics/src/filter_stats.rs``, ``FilterStatsMetrics``).
- ``correct_metrics`` (``fgumi correct --metrics``): ``umi``, ``total_matches``, ``perfect_matches``, ``one_mismatch_matches``, ``two_mismatch_matches``, ``other_matches``, ``fraction_of_matches``, ``representation``; the UMI of only ``N`` counts the unmatched UMIs (``crates/fgumi-metrics/src/correct.rs``, ``UmiCorrectionMetrics``).
- ``clip_metrics`` (``fgumi clip --metrics``): one row per ``read_type`` (``Fragment``, ``ReadOne``, ``ReadTwo``, ``Pair``, ``All``) with ``reads``, ``reads_unmapped``, ``reads_clipped_pre``, ``reads_clipped_post``, ``reads_clipped_five_prime``, ``reads_clipped_three_prime``, ``reads_clipped_overlapping``, ``reads_clipped_extending``, ``bases`` and the same ``bases_clipped_*`` columns (``crates/fgumi-metrics/src/clip.rs``, ``ClippingMetrics``).
"""

from __future__ import annotations

from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import (
    Chart,
    Column,
    Summary,
    Table,
    bars,
    concat,
    from_rows,
    metric_strips,
    number,
    ratio,
    samtools,
    thin,
    ticks,
)

#: The kinds of files that :mod:`.samtools` parses
SAMTOOLS_KINDS = ("bamstats", "flagstats", "idxstats")

#: Family size up to which the family size charts go
MAX_FAMILY_SIZE = 30

#: AB size up to which the duplex family heatmap goes
MAX_DUPLEX_SIZE = 12

#: UMIs per library that the UMI counts table keeps, the most observed; random UMIs have many
MAX_UMIS = 1000

#: The families of each consensus caller, with their labels
FAMILIES = {
    "duplex": {"cs": "CS families", "ss": "SS families", "ds": "DS families"},
    "simplex": {"cs": "CS families", "ss": "SS families"},
}

#: The UMI columns shown in the summary table after the samtools columns;
#: the other UMI counts follow in the row details
UMI_COLUMNS = {
    "accepted_sam_records": bars(),
    "discarded_fraction": ticks(),
    "mean_family_size": ticks(precision=2),
    "frac_raw_reads_used": ticks(),
    "consensus_reads_emitted": bars(),
    "pass_rate": ticks(),
}

#: The yield columns of the full data shown in the summary table, by consensus caller
YIELD_COLUMNS = {
    "duplex": {"ds_duplexes": bars(), "ds_fraction_duplexes": ticks()},
    "simplex": {
        "ss_consensus_families": bars(),
        "mean_ss_family_size": ticks(precision=2),
        "ss_singleton_fraction": ticks(),
    },
}

#: The ``fgumi correct`` columns of the summary table
CORRECT_COLUMNS = {
    "umi_corrected_fraction": ticks(precision=4),
    "umi_unmatched_fraction": ticks(precision=4),
}

#: The ``fgumi correct`` counts of corrected UMIs
MISMATCH_COLUMNS = ("one_mismatch_matches", "two_mismatch_matches", "other_matches")

#: The UMI metrics of the strip chart, with their labels
UMI_METRICS = {
    "discarded_fraction": "discarded records",
    "mean_family_size": "mean family size",
    "frac_raw_reads_used": "raw reads in a consensus",
    "pass_rate": "filter pass rate",
    "ds_fraction_duplexes": "duplex fraction",
    "ss_singleton_fraction": "SS singletons",
    "umi_corrected_fraction": "corrected UMIs",
    "umi_unmatched_fraction": "unmatched UMIs",
}

#: The column of the yield curve and its title, by consensus caller
YIELD_CURVE = {
    "duplex": ("ds_duplexes", "duplexes"),
    "simplex": ("ss_consensus_families", "SS families of consensus size"),
}


def read_metrics(path: str) -> pl.DataFrame:
    """A headered fgumi metrics TSV.
    The types come from all rows, as a fraction column may start with ``0`` or ``1``."""
    return pl.read_csv(path, separator="\t", infer_schema_length=None)


def read_key_values(path: str) -> dict[str, float | int]:
    """The ``key value description`` rows of ``fgumi duplex|simplex --stats`` as key -> number."""
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    return {key: number(value) for key, value in zip(frame["key"], frame["value"])}


def umi_columns(tables: dict[str, pl.DataFrame], stats: dict[str, float | int]) -> dict[str, Any]:
    """The UMI columns of one library's summary row, from its metrics by kind and its consensus stats."""
    grouping = tables["grouping_metrics"].row(0, named=True)
    discarded = sum(value for key, value in grouping.items() if key.startswith("discarded_"))
    groups = tables["group_family_sizes"]
    full = tables["yield_metrics"].sort("fraction").row(-1, named=True)
    row = {
        **grouping,
        "discarded_fraction": ratio(discarded, grouping["accepted_sam_records"] + discarded),
        "mean_family_size": ratio(
            (groups["family_size"] * groups["count"]).sum(), groups["count"].sum()
        ),
        **stats,
        **tables["filter_stats"].row(0, named=True),
        **{key: value for key, value in full.items() if key != "fraction"},
    }
    if "correct_metrics" in tables:
        correct = tables["correct_metrics"]
        total = correct["total_matches"].sum()
        corrected = correct.select(pl.sum_horizontal(MISMATCH_COLUMNS).sum()).item()
        unmatched = correct.filter(pl.col("umi").str.contains("^N+$"))["total_matches"].sum()
        row["umi_corrected_fraction"] = ratio(corrected, total)
        row["umi_unmatched_fraction"] = ratio(unmatched, total)
    return row


def duplex_grid(frame: pl.DataFrame) -> pl.DataFrame:
    """The families by AB and BA size up to :data:`MAX_DUPLEX_SIZE`,
    one row per AB size label (``AB 3``) and BA size.
    The triangle ``ba_size <= ab_size`` is completed with zeros,
    so every library has the same cells,
    and the chart has more AB rows than :data:`MAX_COLORED`, which makes it a heatmap."""
    sizes = [(ab, ba) for ab in range(1, MAX_DUPLEX_SIZE + 1) for ba in range(ab + 1)]
    grid = pl.DataFrame(sizes, schema={"ab_size": pl.Int64, "ba_size": pl.Int64}, orient="row")
    return grid.join(
        frame.select("ab_size", "ba_size", "count"),
        on=["ab_size", "ba_size"],
        how="left",
        maintain_order="left",
    ).select(
        pl.format("AB {}", "ab_size").alias("ab"),
        "ba_size",
        pl.col("count").fill_null(0).alias("families"),
    )


def library_frames(tables: dict[str, pl.DataFrame], consensus: str) -> dict[str, pl.DataFrame]:
    """The rows of one library in each view, by view name."""
    families = FAMILIES[consensus]
    yields = tables["yield_metrics"]
    frames = {
        "family_sizes": tables["family_sizes"]
        .filter(pl.col("family_size") <= MAX_FAMILY_SIZE)
        .select("family_size", *(pl.col(f"{k}_fraction").alias(v) for k, v in families.items()))
        .unpivot(index="family_size", variable_name="kind", value_name="fraction"),
        "group_family_sizes": tables["group_family_sizes"]
        .filter(pl.col("family_size") <= MAX_FAMILY_SIZE)
        .select("family_size", "fraction"),
        "position_group_sizes": tables["position_group_sizes"].select(
            "position_group_size", "fraction"
        ),
        "yield": yields.select("fraction", YIELD_CURVE[consensus][0]),
        "umi_counts": tables["umi_counts"]
        .sort("raw_observations", descending=True, maintain_order=True)
        .head(MAX_UMIS),
    }
    if consensus == "duplex":
        frames["duplex_family_sizes"] = duplex_grid(tables["duplex_family_sizes"])
        frames["duplex_yield"] = yields.select(
            "fraction",
            (pl.col("ds_fraction_duplexes") / pl.col("ds_fraction_duplexes_ideal")).alias(
                "of_ideal"
            ),
        )
    for kind in ("correct_metrics", "clip_metrics"):
        if kind in tables:
            frames[kind] = tables[kind]
    return frames


def _stack(frames: list[pl.DataFrame]) -> pl.DataFrame:
    """The rows of all libraries, the ``library`` column first."""
    return pl.concat(frames, how="diagonal_relaxed").select("library", pl.exclude("library"))


def summarize(files: dict[str, list[str]], labels: list[str], consensus: str) -> Summary:
    bam = samtools.summarize({kind: files[kind] for kind in SAMTOOLS_KINDS}, labels)
    rows: list[dict[str, Any]] = []
    frames: dict[str, list[pl.DataFrame]] = {}
    for index, label in enumerate(labels):
        paths = {kind: files[kind][index] for kind in files if kind not in SAMTOOLS_KINDS}
        stats = read_key_values(paths.pop("consensus_stats"))
        tables = {kind: read_metrics(path) for kind, path in paths.items()}
        rows.append({"library": label, **umi_columns(tables, stats)})
        library = pl.lit(label).alias("library")
        for name, frame in library_frames(tables, consensus).items():
            frames.setdefault(name, []).append(frame.with_columns(library))

    shown = UMI_COLUMNS | YIELD_COLUMNS[consensus]
    if "correct_metrics" in files:
        shown |= CORRECT_COLUMNS
    summary = bam.table.frame.join(
        from_rows(rows, ["library", *shown]), on="library", how="left", maintain_order="left"
    )
    columns = samtools.COLUMNS | {"reads": bars(label="consensus reads")} | shown
    columns |= {
        c: Column(detail=True) for c in summary.columns if c != "library" and c not in columns
    }
    families = FAMILIES[consensus]
    curve, curve_title = YIELD_CURVE[consensus]

    views: list[Table | Chart] = [
        metric_strips(
            summary,
            samtools.METRICS | UMI_METRICS,
            "libraries",
            "Libraries by metric",
            desc="One tick per library on each metric's range: a tick apart from the others "
            "is a library to look at; the tooltip names it.",
        ),
        Chart(
            "family_sizes",
            "Family sizes of the consensus calling",
            concat(frames["family_sizes"], ["library", "kind", "family_size", "fraction"]),
            x="family_size",
            y="fraction",
            facet="kind",
            facet_order=list(families.values()),
            x_title="family size (reads)",
            y_title="fraction of families",
            desc=f"The fraction of the families of each size up to {MAX_FAMILY_SIZE} reads, as "
            f"`fgumi {consensus}` counts them: CS families by coordinate and strand, SS families "
            "by UMI on one strand"
            + (", DS families by UMI on both strands" if "ds" in families else "")
            + ". Most families at size 1 mean that few molecules were read more than once; "
            "families far larger than the consensus needs mean reads beyond saturation.",
        ),
        Chart(
            "group_family_sizes",
            "Family sizes of the UMI grouping",
            concat(frames["group_family_sizes"], ["library", "family_size", "fraction"]),
            x="family_size",
            y="fraction",
            x_title="family size (reads)",
            y_title="fraction of families",
            desc=f"The fraction of the families of `fgumi group` of each size up to "
            f"{MAX_FAMILY_SIZE} reads; the mean family size of the summary counts all sizes.",
        ),
        Chart(
            "position_group_sizes",
            "Families per position",
            thin(
                concat(
                    frames["position_group_sizes"],
                    ["library", "position_group_size", "fraction"],
                ),
                "position_group_size",
                "fraction",
            ),
            x="position_group_size",
            y="fraction",
            x_title="UMI families at the same position",
            y_title="fraction of positions",
            desc="The fraction of the positions (the same start and end coordinates) with each "
            "number of UMI families (`fgumi group`).",
        ),
    ]
    if consensus == "duplex":
        views.append(
            Chart(
                "duplex_family_sizes",
                "Duplex families by AB and BA size",
                concat(frames["duplex_family_sizes"], ["library", "ab", "ba_size", "families"]),
                x="ba_size",
                y="families",
                series="ab",
                facet="library",
                log_y=True,
                x_title="BA size (reads)",
                y_title="families",
                link=None,
                desc="One panel per library: the families by the reads of their larger (AB, rows) "
                f"and smaller (BA, columns) strand, up to {MAX_DUPLEX_SIZE} reads, on a log "
                "color scale. The column of BA size 0 holds the families seen on one strand "
                "only; a duplex needs reads of both strands.",
            )
        )
    views.append(
        Chart(
            "yield",
            "Yield by sampled fraction",
            concat(frames["yield"], ["library", "fraction", curve]),
            x="fraction",
            y=curve,
            x_title="sampled fraction of the read pairs",
            y_title=curve_title,
            desc=f"The {curve_title} of a random sample of the read pairs, by its fraction "
            "(fgumi's yield metrics): a curve that still rises at fraction 1 gains from deeper "
            "sequencing; a flat one is saturated.",
        )
    )
    if consensus == "duplex":
        views.append(
            Chart(
                "duplex_yield",
                "Duplex fraction of the ideal by sampled fraction",
                concat(frames["duplex_yield"], ["library", "fraction", "of_ideal"]),
                x="fraction",
                y="of_ideal",
                x_title="sampled fraction of the read pairs",
                y_title="duplex fraction / ideal",
                desc="The fraction of the DS families that are duplexes over the fraction an "
                "ideal library would have at the same depth (`ds_fraction_duplexes` over "
                "`ds_fraction_duplexes_ideal`): well below 1 means that the two strands of the "
                "molecules were not both captured.",
            )
        )
    views.append(
        Table(
            "umi_counts",
            "UMI counts",
            _stack(frames["umi_counts"]),
            desc=f"The UMIs of the consensus calling, the {MAX_UMIS} most observed of each "
            "library: the raw observations (reads) and the unique observations (families).",
            columns={
                "raw_observations": bars(),
                "unique_observations": bars(),
                "fraction_raw_observations": ticks(precision=4),
                "fraction_unique_observations": ticks(precision=4),
            },
            pinned=("library", "umi"),
        )
    )
    if "correct_metrics" in frames:
        views.append(
            Table(
                "umi_correction",
                "UMI correction",
                _stack(frames["correct_metrics"]),
                desc="`fgumi correct`: the reads whose UMI matched each expected UMI, by "
                "mismatches; the UMI of only `N` counts the reads that matched none.",
                columns={
                    "total_matches": bars(),
                    "fraction_of_matches": ticks(precision=4),
                    "representation": ticks(precision=2),
                },
                pinned=("library", "umi"),
            )
        )
    if "clip_metrics" in frames:
        views.append(
            Table(
                "clipping",
                "Clipping",
                _stack(frames["clip_metrics"]),
                desc="`fgumi clip`: the reads and aligned bases by read type, and those clipped "
                "for each reason.",
                columns={"reads": bars(), "bases": bars()},
                pinned=("library", "read_type"),
            )
        )
    views += [view for view in bam.views if view.name != "libraries"]

    return Summary(
        table=Table(
            "summary",
            "Consensus BAM and UMI QC",
            summary,
            desc="One row per library. "
            "`samtools stats` of the consensus BAM: the consensus reads and their fractions. "
            "`fgumi group`: the raw records accepted for grouping, the fraction discarded and "
            f"the mean family size. `fgumi {consensus}`: the fraction of the raw reads in a "
            "consensus and the consensus reads emitted. `fgumi filter`: the pass rate. "
            "The yield columns are those of the full data. "
            + (
                "`fgumi correct`: the fractions of corrected and of unmatched UMIs. "
                if "correct_metrics" in files
                else ""
            )
            + "The details of a row hold the counts and the reasons of rejection.",
            columns=columns,
        ),
        views=views,
    )
