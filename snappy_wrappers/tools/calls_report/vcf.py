# -*- coding: utf-8 -*-
"""Reading a VCF as text: the header's fields and one record per ALT allele.

The values stay the text of the VCF; the tiers of the calls report parse what they show.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field

_HEADER_LINE = re.compile(r"##(INFO|FORMAT|contig)=<(.*)>$")
_ATTRIBUTE = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,]*)')


@dataclass
class Record:
    """One ALT allele of one VCF record."""

    chromosome: str
    position: int
    id: str
    ref: str
    alt: str
    #: Index of the allele among the record's ALT alleles, for ``AD`` and the annotations
    alt_index: int
    #: Number of ALT alleles of the record
    alts: int
    quality: str
    filter: str
    #: INFO field -> text (``""`` for a flag)
    info: dict[str, str]
    #: Sample -> FORMAT field -> text
    samples: dict[str, dict[str, str]]


@dataclass
class Vcf:
    """The header fields and the records of a VCF."""

    samples: list[str]
    #: INFO ID -> header attributes (``Number``, ``Type``, ``Description``, ...), in header order
    info: dict[str, dict[str, str]] = field(default_factory=dict)
    #: FORMAT ID -> header attributes, in header order
    formats: dict[str, dict[str, str]] = field(default_factory=dict)
    #: Contig names in header order
    contigs: list[str] = field(default_factory=list)
    records: list[Record] = field(default_factory=list)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def header_attributes(text: str) -> dict[str, str]:
    """The ``key=value`` attributes of a header line's ``<...>`` part, quotes removed."""
    attributes = {}
    for key, value in _ATTRIBUTE.findall(text):
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1].replace('\\"', '"')
        attributes[key] = value
    return attributes


def read_vcf(path: str) -> Vcf:
    """Read ``path`` into a :class:`Vcf`, one :class:`Record` per ALT allele."""
    vcf = Vcf(samples=[])
    with _open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                if match := _HEADER_LINE.match(line):
                    kind, text = match.groups()
                    attributes = header_attributes(text)
                    if kind == "INFO":
                        vcf.info[attributes["ID"]] = attributes
                    elif kind == "FORMAT":
                        vcf.formats[attributes["ID"]] = attributes
                    else:
                        vcf.contigs.append(attributes["ID"])
                continue
            if line.startswith("#"):
                vcf.samples = line.split("\t")[9:]
                continue
            fields = line.split("\t")
            chromosome, position, id_, ref, alt, quality, filter_, info_text = fields[:8]
            info: dict[str, str] = {}
            for entry in info_text.split(";"):
                if entry and entry != ".":
                    key, _, value = entry.partition("=")
                    info[key] = value
            samples: dict[str, dict[str, str]] = {}
            if len(fields) > 9:
                keys = fields[8].split(":")
                for sample, text in zip(vcf.samples, fields[9:]):
                    samples[sample] = dict(zip(keys, text.split(":")))
            alts = alt.split(",")
            for index, allele in enumerate(alts):
                vcf.records.append(
                    Record(
                        chromosome=chromosome,
                        position=int(position),
                        id=id_,
                        ref=ref,
                        alt=allele,
                        alt_index=index,
                        alts=len(alts),
                        quality=quality,
                        filter=filter_,
                        info=info,
                        samples=samples,
                    )
                )
    return vcf
