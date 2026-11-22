# -*- coding: utf-8 -*-
"""The calls of one VCF as a datavzrd site for the report, and the same table as a result file.

One row per call (per ALT allele); the views of the site are in :mod:`.views`:
one per event of a varlociraptor scenario, ranked by posterior, rarity and impact,
and all calls in genomic order.
The columns follow what the VCF carries:
with varlociraptor's ``PROB_*`` fields the posterior probability of every event and per sample the observations (:mod:`.varlociraptor`),
otherwise the genotype, allele frequency, alt reads,
depth and genotype quality of each sample (:mod:`.samples`);
with an ``ANN`` or ``CSQ`` annotation the gene, consequence, impact and HGVS of the picked entry,
and with mehari's frequency and ClinVar fields the gnomAD frequency and the ClinVar classification (:mod:`.annotation`).
A field the VCF lacks gives no column; nothing is required.
Linkouts (gnomAD, ClinVar, dbSNP, Ensembl, GeneCards, OMIM,
varlociraptor-inspect) are datavzrd links on the position, gene, ID, ClinVar and inspect columns;
the build-dependent ones need the assembly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import polars as pl

from snappy_wrappers.tools.qc_summary import Column, Summary, ticks

from . import annotation, pedigree, samples, varlociraptor, views
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
    #: The members of the entity's pedigree (see :func:`pedigree.members_of`); empty for none
    pedigree: tuple[dict[str, Any], ...] = ()


def variant_type(ref: str, alt: str, info: dict[str, str]) -> str:
    if "SVTYPE" in info:
        return info["SVTYPE"]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INS" if len(alt) > len(ref) else "DEL"


def sv_end(record: Record) -> int:
    """The last reference base of the variant: ``END`` of a structural variant, else by the REF length."""
    end = record.info.get("END")
    return int(end) if end else record.position + len(record.ref) - 1


def sv_length(record: Record) -> int | None:
    """The length of a structural variant: ``SVLEN`` (of the allele), else the span to ``END``;
    ``None`` for a breakend, which has neither."""
    values = record.info.get("SVLEN", "").split(",")
    value = values[record.alt_index] if len(values) > 1 else values[0]
    if value not in ("", "."):
        return abs(int(value))
    end = record.info.get("END")
    return int(end) - record.position if end else None


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
    """The rows of the calls table and the rendering of its columns, in column order:
    the core (with the length and the end of structural variants), the event probabilities,
    the annotation, the frequencies and ClinVar, the samples,
    then ID, quality, filter, hints and the inspect link."""
    events = varlociraptor.events(vcf)
    annotated = annotation.annotation_of(vcf)
    frequencies = annotation.frequency_columns(vcf, annotated)
    clinvar = annotation.has_clinvar(vcf, annotated)
    links = annotation.linkouts(set(options.linkouts), options.assembly)
    inspect = events and "varlociraptor_inspect" in options.linkouts
    members = pedigree.members_of(list(options.pedigree), vcf.samples)
    with_pedigree = bool(events) and len(members) > 1 and any(m.sample for m in members)
    # By the records: varlociraptor declares the SV fields in every header
    structural = any("SVTYPE" in record.info for record in vcf.records)
    rows = []
    for record in vcf.records:
        row = _core(record)
        if structural:
            row["SV length"] = sv_length(record)
            row["end"] = sv_end(record)
        if links["position"]:
            row["gnomad_chromosome"] = annotation.gnomad_chromosome(record.chromosome)
        if events:
            row["event"] = ""
            row |= varlociraptor.probabilities(record, events)
            row["event"] = views.event_of(row, events)
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
    names = list(CORE)
    if structural:  # the length after the type; the end, hidden, for the IGV views
        names[3:3] = ["SV length", "end"]
    columns: dict[str, Column] = {name: Column() for name in names}
    if structural:
        columns["SV length"] = Column(precision=0)
        columns["end"] = Column(display="hidden")
    columns["position"] = Column(precision=0, links=links["position"] or None)
    if links["position"]:
        columns["gnomad_chromosome"] = Column(display="hidden")
    if events:
        columns["event"] = views.event_column(events)
        columns |= varlociraptor.probability_columns(vcf)
    if with_pedigree:  # first in the row details
        for row in rows:
            row["pedigree"] = ""
        columns["pedigree"] = pedigree.column(members)
    if annotated:
        columns |= annotation.columns(annotated, links["gene"])
    columns |= annotation.frequency_column_specs(frequencies)
    if clinvar:
        columns |= annotation.clinvar_columns(links["ClinVar"])
    if events:
        columns |= varlociraptor.sample_columns(vcf, rows, on_demand=with_pedigree)
    else:
        columns |= samples.columns(vcf)
    columns["id"] = Column(display="available", links=links["id"] or None)
    columns["quality"] = ticks(precision=1, display="available")
    filters = {row["filter"] for row in rows} - {"PASS", "."}
    columns["filter"] = Column(display=None if filters else "available")
    if "HINTS" in vcf.info:
        columns["hints"] = Column(display="available")
    if inspect:
        columns |= varlociraptor.inspect_columns()
    for name in ("other annotations", "hints"):  # last among the columns
        if name in columns:
            columns[name] = columns.pop(name)
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


@dataclass
class Description:
    """A view description for datavzrd, which renders it as Markdown:
    a lead sentence, then titled bullet lists (``Columns``, ``Links``, ``Navigation``)."""

    lead: str
    sections: dict[str, list[str]] = field(default_factory=dict)

    def with_lead(self, lead: str, **extra: list[str]) -> "Description":
        """A copy with another lead and further bullets per section."""
        sections = {title: list(items) for title, items in self.sections.items()}
        for title, items in extra.items():
            sections.setdefault(title.capitalize(), []).extend(items)
        return Description(lead, sections)

    def markdown(self) -> str:
        parts = [self.lead]
        for title, items in self.sections.items():
            if items:
                parts.append(f"**{title}**\n" + "\n".join(f"- {item}" for item in items))
        return "\n\n".join(parts)


def describe(vcf: Vcf, columns: dict[str, Column], options: Options) -> Description:
    """What the columns of the table show, as the bullets every view shares."""
    column_notes: list[str] = []
    if varlociraptor.events(vcf):
        column_notes += varlociraptor.COLUMN_NOTES
    if "pedigree" in columns:
        column_notes.insert(
            1,
            "`pedigree` (row details, first): the pedigree with, per member, the call's allele "
            "frequency as a bar with its 90 % credible interval, below it the reads (for the "
            "alternative allele to the right in blue, for the reference to the left in gray, "
            "darker for stronger evidence, pale for uncertain mapping, scaled to the member "
            "with the most reads, whose width the gray background shows), the depth and the "
            "genotype. "
            "Squares in aqua are male, circles in orange female, diamonds in gray of unknown "
            "sex; an affected member is filled in its color, an unaffected one white, one of "
            "unknown status light gray; a dashed outline has no data; a ring marks the member "
            "of a per-member event such as `denovo_<member>` (the legend under the pedigree). "
            "The per-sample `AF`, `DP`, `GT`, observations and density columns are on demand.",
        )
    else:
        column_notes.append(
            "Per sample: genotype, allele frequency, reads for the alternative allele, depth and "
            "genotype quality, as far as the VCF has them."
        )
    if "gene" in columns:
        column_notes.append(
            "`gene`, `consequence`, `impact`, `hgvs.c`, `hgvs.p`: the picked annotation, "
            "the first of the highest impact for the allele; transcript, biotype, exon and the "
            "other annotations on demand."
        )
    if "gnomAD AF" in columns:
        column_notes.append(
            "`gnomAD AF`: the larger of the exomes and genomes frequencies, empty when not "
            "observed; the single frequencies and the homozygote count on demand."
        )
    if "ClinVar" in columns:
        column_notes.append("`ClinVar`: the germline classifications; the VCV accession on demand.")
    if "SV length" in columns:
        column_notes.append(
            "`SV length`: the length of a structural variant, `SVLEN` or else the span to `END`; "
            "a breakend has none."
        )
    on_demand = [name for name in ("id", "quality", "filter") if name in columns]
    if on_demand:
        column_notes.append(", ".join(f"`{name}`" for name in on_demand) + ": on demand.")
    if "hints" in columns:
        column_notes.append("`hints`: varlociraptor's hints, on demand.")
    links: list[str] = []
    for name in ("position", "gene", "id", "ClinVar"):
        if name in columns and columns[name].links:
            links.append(f"`{name}`: " + ", ".join(columns[name].links))
    if "inspect" in columns:
        links.append(
            "`inspect`: the observations of the call in varlociraptor-inspect, "
            "the samples numbered in the order of the VCF (the order of the sample columns)."
        )
    if options.linkouts and options.assembly not in annotation.ASSEMBLIES:
        links.append(
            "The variants carry no genome build (no annotation task in the chain), "
            "so the gnomAD, ClinVar, dbSNP and Ensembl links by position and gene are left out."
        )
    return Description("", {"Columns": column_notes, "Links": links})


def summarize(path: str, options: Options = Options()) -> Summary:
    """The calls of one VCF as the views of a site, ready to render."""
    vcf = read_vcf(path)
    rows, columns = rows_and_columns(vcf, options)
    return views.build(
        vcf,
        frame_of(rows, columns),
        columns,
        options.page_size,
        options.max_in_memory_rows,
        describe(vcf, columns, options),
    )
