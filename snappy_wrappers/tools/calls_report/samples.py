# -*- coding: utf-8 -*-
"""The per-sample columns of a VCF without varlociraptor's fields: genotype, allele frequency, alt reads, depth and genotype quality from the FORMAT fields the header declares."""

from __future__ import annotations

from typing import Any

from snappy_wrappers.tools.qc_summary import Column, bars, ticks

from .vcf import Record, Vcf


def _number(text: str | None) -> float | None:
    if text in (None, "", "."):
        return None
    return float(text)


def _allele_depth(values: dict[str, str], alt_index: int) -> tuple[int | None, int | None]:
    """The reads of the ALT allele and of all alleles from ``AD``; ``None`` without it."""
    counts = [_number(c) for c in values.get("AD", "").split(",")]
    if len(counts) <= alt_index + 1 or any(c is None for c in counts):
        return None, None
    return int(counts[alt_index + 1]), int(sum(counts))


def values(record: Record, vcf: Vcf) -> dict[str, Any]:
    """The per-sample values of ``record``, only for the FORMAT fields the header declares."""
    row: dict[str, Any] = {}
    for sample, fields in record.samples.items():
        if "GT" in vcf.formats:
            row[f"{sample}: GT"] = fields.get("GT", "")
        alt_reads, depth = (
            _allele_depth(fields, record.alt_index) if "AD" in vcf.formats else (None, None)
        )
        if "AF" in vcf.formats:
            frequencies = fields.get("AF", "").split(",")
            text = (
                frequencies[record.alt_index]
                if len(frequencies) > record.alt_index
                else frequencies[0]
            )
            row[f"{sample}: AF"] = _number(text)
        elif "AD" in vcf.formats:
            row[f"{sample}: AF"] = alt_reads / depth if depth else None
        if "AD" in vcf.formats:
            row[f"{sample}: alt reads"] = alt_reads
        if "DP" in vcf.formats:
            dp = _number(fields.get("DP"))
            row[f"{sample}: DP"] = None if dp is None else int(dp)
        if "GQ" in vcf.formats:
            row[f"{sample}: GQ"] = _number(fields.get("GQ"))
    return row


def columns(vcf: Vcf) -> dict[str, Column]:
    """The rendering of the per-sample columns, in order."""
    result: dict[str, Column] = {}
    for sample in vcf.samples:
        if "GT" in vcf.formats:
            result[f"{sample}: GT"] = Column()
        if "AF" in vcf.formats or "AD" in vcf.formats:
            result[f"{sample}: AF"] = ticks(precision=3, domain=(0, 1))
        if "AD" in vcf.formats:
            result[f"{sample}: alt reads"] = bars()
        if "DP" in vcf.formats:
            result[f"{sample}: DP"] = bars()
        if "GQ" in vcf.formats:
            result[f"{sample}: GQ"] = ticks(precision=0, display="available")
    return result
