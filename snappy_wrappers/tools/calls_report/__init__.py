# -*- coding: utf-8 -*-
"""The calls of one VCF as a datavzrd site for the report, and the same table as a result file.

One row per call (per ALT allele).
The columns follow what the VCF carries: with varlociraptor's ``PROB_*`` fields the posterior probability of every event and per sample the observations (:mod:`.varlociraptor`), otherwise the genotype, allele frequency, alt reads, depth and genotype quality of each sample (:mod:`.samples`); with an ``ANN`` or ``CSQ`` annotation the gene, consequence, impact and HGVS of the picked entry, and with mehari's frequency and ClinVar fields the gnomAD frequency and the ClinVar classification (:mod:`.annotation`).
A field the VCF lacks gives no column; nothing is required.
Linkouts (gnomAD, ClinVar, dbSNP, Ensembl, GeneCards, OMIM, varlociraptor-inspect) are datavzrd links on the position, gene, ID, ClinVar and inspect columns; the build-dependent ones need the assembly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import Column, Summary, Table, ticks

from . import annotation, samples, varlociraptor
from .vcf import Record, Vcf, read_vcf

#: The first columns of every table
CORE = ("chromosome", "position", "type", "ref", "alt")


@dataclass(frozen=True)
class Options:
    """The settings of the report from the task's ``calls_report`` section."""

    page_size: int = 25
    max_in_memory_rows: int | None = None
    #: The linkouts to render (names of ``snappy_pipeline.models.calls_report.Linkout``)
    linkouts: tuple[str, ...] = ()
    #: The genome build of the variants (``GRCh37``, ``GRCh38``), for the build-dependent linkouts
    assembly: str | None = None


def variant_type(ref: str, alt: str, info: dict[str, str]) -> str:
    if "SVTYPE" in info:
        return info["SVTYPE"]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INS" if len(alt) > len(ref) else "DEL"


def _core(record: Record) -> dict[str, Any]:
    return {
        "chromosome": record.chromosome,
        "position": record.position,
        "type": variant_type(record.ref, record.alt, record.info),
        "ref": record.ref,
        "alt": record.alt,
    }


def _trailing(record: Record, vcf: Vcf) -> dict[str, Any]:
    row: dict[str, Any] = {
        "id": "" if record.id == "." else record.id,
        "quality": None if record.quality == "." else float(record.quality),
        "filter": record.filter,
    }
    if "HINTS" in vcf.info:
        row["hints"] = record.info.get("HINTS", "")
    return row


def rows_and_columns(
    vcf: Vcf, options: Options = Options()
) -> tuple[list[dict[str, Any]], dict[str, Column]]:
    """The rows of the calls table and the rendering of its columns, in column order: the
    core, the event probabilities, the annotation, the frequencies and ClinVar, the samples,
    then ID, quality, filter, hints and the inspect link."""
    events = varlociraptor.events(vcf)
    annotated = annotation.annotation_of(vcf)
    frequencies = annotation.frequency_columns(vcf, annotated)
    clinvar = annotation.has_clinvar(vcf, annotated)
    links = annotation.linkouts(set(options.linkouts), options.assembly)
    inspect = events and "varlociraptor_inspect" in options.linkouts
    rows = []
    for record in vcf.records:
        row = _core(record)
        if links["position"]:
            row["gnomad_chromosome"] = annotation.gnomad_chromosome(record.chromosome)
        if events:
            row |= varlociraptor.probabilities(record, events)
        if annotated:
            row |= annotation.values(record, annotated)
        row |= annotation.frequency_values(record, frequencies, annotated)
        if clinvar:
            row |= annotation.clinvar_values(record, annotated)
        if events:
            row |= varlociraptor.sample_values(record, vcf.samples)
        else:
            row |= samples.values(record, vcf)
        row |= _trailing(record, vcf)
        if inspect:
            row |= varlociraptor.inspect_values(record, events, vcf.samples)
        rows.append(row)
    columns: dict[str, Column] = {name: Column() for name in CORE}
    columns["position"] = Column(precision=0, links=links["position"] or None)
    if links["position"]:
        columns["gnomad_chromosome"] = Column(display="hidden")
    if events:
        columns |= varlociraptor.probability_columns(vcf)
    if annotated:
        columns |= annotation.columns(annotated, links["gene"])
    columns |= annotation.frequency_column_specs(frequencies)
    if clinvar:
        columns |= annotation.clinvar_columns(links["ClinVar"])
    if events:
        columns |= varlociraptor.sample_columns(vcf, rows)
    else:
        columns |= samples.columns(vcf)
    columns["id"] = Column(display="available", links=links["id"] or None)
    columns["quality"] = ticks(precision=1, display="available")
    filters = {row["filter"] for row in rows} - {"PASS", "."}
    columns["filter"] = Column(display=None if filters else "available")
    if "HINTS" in vcf.info:
        columns["hints"] = Column(display="detail")
    if inspect:
        columns |= varlociraptor.inspect_columns()
    for name in ("id", "quality", "filter"):  # a column the VCF does not fill is left out
        if rows and all(row[name] in (None, "", ".") for row in rows):
            del columns[name]
            for row in rows:
                del row[name]
    return rows, columns


def frame_of(rows: list[dict[str, Any]], columns: dict[str, Column]) -> pl.DataFrame:
    """A frame of ``rows`` with the columns in the order of ``columns``."""
    if rows:
        frame = pl.from_dicts(rows, infer_schema_length=None)
    else:
        frame = pl.DataFrame(schema={column: pl.Null for column in columns})
    return frame.select(*(c for c in columns if c in frame.columns))


def summarize(path: str, options: Options = Options()) -> Summary:
    """The calls table of one VCF, ready to render."""
    vcf = read_vcf(path)
    rows, columns = rows_and_columns(vcf, options)
    desc = "One row per call. "
    if varlociraptor.events(vcf):
        desc += varlociraptor.DESCRIPTION
    else:
        desc += (
            "Per sample: the genotype, the allele frequency, the reads for the alternative "
            "allele, the depth and the genotype quality, as far as the VCF has them. "
        )
    if "gene" in columns:
        desc += (
            "Gene, consequence, impact and HGVS are those of the picked annotation (the first "
            "of the highest impact for the allele); the other annotations are in the details. "
        )
    if "gnomAD AF" in columns:
        desc += "gnomAD AF is the larger of the exomes and genomes frequencies; empty means not observed. "
    if options.linkouts and options.assembly not in annotation.ASSEMBLIES:
        desc += (
            "The variants carry no genome build (no annotation task in the chain), so the "
            "gnomAD, ClinVar, dbSNP and Ensembl linkouts by position and gene are left out. "
        )
    if "inspect" in columns:
        desc += (
            "The inspect link opens the observations of the call in varlociraptor-inspect, "
            "the samples numbered in the order of the VCF (the order of the sample columns). "
        )
    desc += "The ID, quality and filter columns are available on demand."
    table = Table(
        "all",
        "all calls",
        frame_of(rows, columns),
        desc=desc,
        columns=columns,
        pinned=("chromosome", "position"),
        page_size=options.page_size,
        max_in_memory_rows=options.max_in_memory_rows,
    )
    return Summary(table=table, views=[])
