# -*- coding: utf-8 -*-
"""The annotation columns of a VCF and the linkouts.

The consequence of the picked annotation comes from ``ANN`` (mehari, Jannovar, SnpEff) or ``CSQ`` (VEP): the entries for the row's allele, the first of the highest impact.
Population frequencies come from the fields mehari adds with ``--frequencies`` (``gnomad_exomes_*``, ``gnomad_genomes_*``, which hold gnomAD mtDNA on chrM, and ``helix_*``: allele numbers and carrier counts, from which the allele frequency is computed) or from VEP's ``gnomADe_AF`` and ``gnomADg_AF``; the ClinVar classification from mehari's ``clinvar_germline_classification`` and ``clinvar_vcv`` or VEP's ``CLIN_SIG``.
A VCF without a field gets no column.
The linkouts are datavzrd ``link-to-url`` entries; the ones that depend on the genome build (gnomAD, ClinVar and dbSNP by position, Ensembl) need the assembly, which the annotation task tags its variants with.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any
from urllib.parse import unquote

from snappy_wrappers.tools.qc_summary import INK, MUTED, STATUS, Column, pills, ticks

from .vcf import Record, Vcf

#: The impacts, highest first, and their colors (status colors; the label is always shown)
IMPACTS = ("HIGH", "MODERATE", "LOW", "MODIFIER")
IMPACT_COLORS = (STATUS["critical"], STATUS["serious"], STATUS["warning"], MUTED)

#: ClinVar's germline classifications and their colors: pathogenic on the status colors,
#: benign in blue, conflicting in violet, so that no red/green pair carries the meaning
CLINVAR = (
    ("Pathogenic", STATUS["critical"]),
    ("Pathogenic/Likely pathogenic", STATUS["critical"]),
    ("Likely pathogenic", STATUS["serious"]),
    ("Uncertain significance", STATUS["warning"]),
    ("Likely benign", "#8fb8ea"),
    ("Benign/Likely benign", INK),
    ("Benign", INK),
    ("Conflicting classifications of pathogenicity", "#6250d6"),
)

#: Genome builds the linkouts know
ASSEMBLIES = ("GRCh37", "GRCh38")

#: Table column -> (``ANN`` column, ``CSQ`` column)
FIELDS = {
    "gene": ("Gene_Name", "SYMBOL"),
    "gene_id": ("Gene_ID", "Gene"),
    "hgnc_id": ("HGNC_ID", "HGNC_ID"),
    "consequence": ("Annotation", "Consequence"),
    "impact": ("Annotation_Impact", "IMPACT"),
    "transcript": ("Feature_ID", "Feature"),
    "biotype": ("Transcript_BioType", "BIOTYPE"),
    "exon": ("Rank", "EXON"),
    "hgvs.c": ("HGVS.c", "HGVSc"),
    "hgvs.p": ("HGVS.p", "HGVSp"),
}

#: The frequency fields of mehari, by table column: (allele number, hom, het, hemi)
MEHARI_FREQUENCIES = {
    "gnomAD exomes AF": (
        "gnomad_exomes_an",
        "gnomad_exomes_hom",
        "gnomad_exomes_het",
        "gnomad_exomes_hemi",
    ),
    "gnomAD genomes AF": (
        "gnomad_genomes_an",
        "gnomad_genomes_hom",
        "gnomad_genomes_het",
        "gnomad_genomes_hemi",
    ),
    "HelixMtDb AF": ("helix_an", "helix_hom", "helix_het", None),
}

#: The frequency columns of VEP, by table column
VEP_FREQUENCIES = {"gnomAD exomes AF": "gnomADe_AF", "gnomAD genomes AF": "gnomADg_AF"}

INSPECT_URL = "https://varlociraptor.github.io/varlociraptor-inspect/"


@dataclass(frozen=True)
class Annotation:
    """The annotation field of the header: its key and the columns of an entry."""

    key: str
    columns: tuple[str, ...]

    @property
    def allele(self) -> int | None:
        return self.columns.index("Allele") if "Allele" in self.columns else None

    def fields(self) -> dict[str, int]:
        """Table column -> index of its value in an entry, for the columns the header has."""
        source = 0 if self.key == "ANN" else 1
        result = {}
        for column, names in FIELDS.items():
            if names[source] in self.columns:
                result[column] = self.columns.index(names[source])
        return result

    def entries(self, record: Record) -> list[list[str]]:
        """The entries of ``record`` for its allele, percent-decoded; all entries when none
        names the allele (VEP minimises alleles)."""
        text = record.info.get(self.key, "")
        entries = [[unquote(v) for v in entry.split("|")] for entry in text.split(",") if entry]
        if self.allele is not None:
            matching = [e for e in entries if len(e) > self.allele and e[self.allele] == record.alt]
            if matching:
                entries = matching
        return entries


def annotation_of(vcf: Vcf) -> Annotation | None:
    """The ``ANN`` or ``CSQ`` annotation the header declares, else ``None``."""
    patterns = {"ANN": r"Functional annotations: ?'(.*)'", "CSQ": r"Format: (.*)"}
    for key, pattern in patterns.items():
        description = vcf.info.get(key, {}).get("Description", "")
        if match := re.search(pattern, description):
            return Annotation(key, tuple(c.strip() for c in match.group(1).split("|")))
    return None


def impact_rank(impact: str | None) -> int:
    return IMPACTS.index(impact) if impact in IMPACTS else len(IMPACTS)


def pick(entries: list[list[str]], impact_index: int | None) -> list[str] | None:
    """The entry shown: the first of the highest impact."""
    if not entries:
        return None
    if impact_index is None:
        return entries[0]
    return min(
        entries,
        key=lambda e: impact_rank(e[impact_index] if len(e) > impact_index else None),
    )


def values(record: Record, annotation: Annotation) -> dict[str, Any]:
    """The annotation values of ``record``: the picked entry's fields, the helper ids for the
    links and the other entries as ``transcript: consequence``."""
    fields = annotation.fields()
    entries = annotation.entries(record)
    chosen = pick(entries, fields.get("impact"))
    row: dict[str, Any] = {}
    for column, index in fields.items():
        value = chosen[index] if chosen and len(chosen) > index else ""
        if column == "gene_id":
            value = value.split(".")[0]
        elif column == "hgnc_id":
            value = value.removeprefix("HGNC:")
        row[column] = value
    if "transcript" in fields and "consequence" in fields:
        others = [e for e in entries if e is not chosen]
        row["other annotations"] = "; ".join(
            f"{e[fields['transcript']]}: {e[fields['consequence']]}"
            for e in others
            if len(e) > max(fields["transcript"], fields["consequence"])
        )
    return row


def columns(annotation: Annotation, links: dict[str, str]) -> dict[str, Column]:
    """The rendering of the annotation columns, in order; ``links`` are the gene linkouts."""
    fields = annotation.fields()
    result: dict[str, Column] = {}
    if "gene" in fields:
        result["gene"] = Column(links=links or None)
    if "consequence" in fields:
        result["consequence"] = Column(ellipsis=32)
    if "impact" in fields:
        result["impact"] = pills(IMPACTS, IMPACT_COLORS)
    for name in ("hgvs.c", "hgvs.p"):
        if name in fields:
            result[name] = Column(ellipsis=28)
    for name in ("transcript", "biotype", "exon"):
        if name in fields:
            result[name] = Column(display="available")
    if "other annotations" in result or ("transcript" in fields and "consequence" in fields):
        result["other annotations"] = Column(display="detail")
    for name in ("gene_id", "hgnc_id"):
        if name in fields:
            result[name] = Column(display="hidden")
    return result


def _integer(text: str | None) -> int | None:
    return None if text in (None, "", ".") else int(text)


def _frequency(info: dict[str, str], names: tuple[str | None, ...]) -> float | None:
    an_name, hom_name, het_name, hemi_name = names
    an = _integer(info.get(an_name))
    if an is None:
        return None
    if an == 0:
        return 0.0
    hom, het = _integer(info.get(hom_name)) or 0, _integer(info.get(het_name)) or 0
    hemi = (_integer(info.get(hemi_name)) or 0) if hemi_name else 0
    return (2 * hom + het + hemi) / an


def frequency_columns(vcf: Vcf, annotation: Annotation | None) -> dict[str, str | tuple]:
    """Table column -> source of a frequency the VCF carries: mehari's fields, or VEP's column."""
    sources: dict[str, str | tuple] = {}
    for column, names in MEHARI_FREQUENCIES.items():
        if names[0] in vcf.info:
            sources[column] = names
    if annotation and annotation.key == "CSQ":
        for column, name in VEP_FREQUENCIES.items():
            if name in annotation.columns:
                sources[column] = name
    return sources


def frequency_values(
    record: Record, sources: dict[str, str | tuple], annotation: Annotation | None
) -> dict[str, Any]:
    """The frequencies of ``record``, and ``gnomAD AF``, the largest of exomes and genomes."""
    row: dict[str, Any] = {}
    chosen = None
    if annotation and annotation.key == "CSQ":
        chosen = pick(annotation.entries(record), annotation.fields().get("impact"))
    for column, source in sources.items():
        if isinstance(source, tuple):
            row[column] = _frequency(record.info, source)
        else:
            index = annotation.columns.index(source)
            text = chosen[index] if chosen and len(chosen) > index else ""
            row[column] = float(text) if text else None
    observed = [row[c] for c in ("gnomAD exomes AF", "gnomAD genomes AF") if row.get(c) is not None]
    if any(c in sources for c in ("gnomAD exomes AF", "gnomAD genomes AF")):
        row["gnomAD AF"] = max(observed) if observed else None
    if "gnomad_exomes_an" in record.info or "gnomad_genomes_an" in record.info:
        row["gnomAD hom"] = sum(
            _integer(record.info.get(name)) or 0
            for name in ("gnomad_exomes_hom", "gnomad_genomes_hom")
        )
    return row


def frequency_column_specs(sources: dict[str, str | tuple]) -> dict[str, Column]:
    result: dict[str, Column] = {}
    if any(c in sources for c in ("gnomAD exomes AF", "gnomAD genomes AF")):
        result["gnomAD AF"] = ticks(precision=4, domain=(0, 1), scale="symlog")
    for column in sources:
        result[column] = ticks(precision=4, domain=(0, 1), scale="symlog", display="available")
    if any(isinstance(s, tuple) and s[0].startswith("gnomad") for s in sources.values()):
        result["gnomAD hom"] = Column(precision=0, display="available")
    return result


def has_clinvar(vcf: Vcf, annotation: Annotation | None) -> bool:
    if "clinvar_germline_classification" in vcf.info:
        return True
    return bool(annotation and annotation.key == "CSQ" and "CLIN_SIG" in annotation.columns)


def clinvar_values(record: Record, annotation: Annotation | None) -> dict[str, Any]:
    """``ClinVar``: the germline classifications; ``ClinVar VCV``: the accessions; ``clinvar_id``:
    the digits of the first accession, for the link."""
    row: dict[str, Any] = {"ClinVar": "", "ClinVar VCV": "", "clinvar_id": ""}
    if "clinvar_germline_classification" in record.info:
        row["ClinVar"] = record.info["clinvar_germline_classification"].replace("%2C", ",")
        vcvs = record.info.get("clinvar_vcv", "")
        row["ClinVar VCV"] = vcvs
        first = vcvs.split(",")[0]
        row["clinvar_id"] = first.split(".")[0].removeprefix("VCV").lstrip("0") if first else ""
    elif annotation and "CLIN_SIG" in annotation.columns:
        chosen = pick(annotation.entries(record), annotation.fields().get("impact"))
        index = annotation.columns.index("CLIN_SIG")
        text = chosen[index] if chosen and len(chosen) > index else ""
        row["ClinVar"] = ",".join(
            part.replace("_", " ").capitalize() for part in text.split("&") if part
        )
    return row


def clinvar_columns(links: dict[str, str]) -> dict[str, Column]:
    domain, colors = zip(*CLINVAR)
    return {
        "ClinVar": pills(domain, colors, ellipsis=24, links=links or None),
        "ClinVar VCV": Column(display="available"),
        "clinvar_id": Column(display="hidden"),
    }


def gnomad_chromosome(chromosome: str) -> str:
    name = chromosome.removeprefix("chr")
    return "M" if name == "MT" else name


def linkouts(selected: set[str], assembly: str | None) -> dict[str, dict[str, str]]:
    """The linkouts by table column (``position``, ``gene``, ``id``, ``ClinVar``), for the
    ``selected`` names; the build-dependent ones only with a known ``assembly``."""
    build = assembly in ASSEMBLIES
    result: dict[str, dict[str, str]] = {"position": {}, "gene": {}, "id": {}, "ClinVar": {}}
    if "gnomad" in selected and build:
        dataset = "gnomad_r4" if assembly == "GRCh38" else "gnomad_r2_1"
        result["position"]["gnomAD"] = (
            "https://gnomad.broadinstitute.org/variant/{gnomad_chromosome}-{position}-{ref}-{alt}"
            f"?dataset={dataset}"
        )
    if "clinvar" in selected and build:
        field = "chrpos38" if assembly == "GRCh38" else "chrpos37"
        result["position"]["ClinVar"] = (
            "https://www.ncbi.nlm.nih.gov/clinvar/?term={gnomad_chromosome}[chr]+AND+{position}"
            f"[{field}]"
        )
        result["ClinVar"]["ClinVar"] = (
            "https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
        )
    if "dbsnp" in selected:
        result["id"]["dbSNP"] = "https://www.ncbi.nlm.nih.gov/snp/{value}"
        if build:
            field = "Base Position" if assembly == "GRCh38" else "Base Position Previous"
            result["position"]["dbSNP"] = (
                "https://www.ncbi.nlm.nih.gov/snp/?term={gnomad_chromosome}[Chromosome]+AND+"
                f"{{position}}[{field}]"
            )
    if "ensembl" in selected and build:
        host = "www.ensembl.org" if assembly == "GRCh38" else "grch37.ensembl.org"
        result["gene"]["Ensembl"] = f"https://{host}/Homo_sapiens/Gene/Summary?g={{gene_id}}"
    if "genecards" in selected:
        result["gene"]["GeneCards"] = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={value}"
    if "omim" in selected:
        result["gene"]["OMIM"] = "https://www.omim.org/search?search={value}"
    return result
