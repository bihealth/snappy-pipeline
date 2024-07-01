#!/usr/bin/env python3
"""Helper tool for merging exome CNV results."""

import argparse
import contextlib
import logging
import os
import sys
import typing
from statistics import mean

import attr
import logzero
import ncls
import pandas as pd
import vcfpy
from logzero import logger

#: Source program is GATK gCNV
SOURCE_GATK_GCNV = "GCNV"

#: Type of the CNV is DEL.
CNV_DEL = "DEL"
#: Type of the CNV is DUP.
CNV_DUP = "DUP"

#: Mapping from VCF "source" header value to internal representation.
SOURCE_MAP = {
    "PostprocessGermlineCNVCalls": SOURCE_GATK_GCNV,
    "JointGermlineCNVSegmentation": SOURCE_GATK_GCNV,
}


class UnionFind:
    """Union-Find (disjoint set) data structure allowing to address by vertex name"""

    def __init__(self, vertex_names):
        #: Node name to id mapping
        self._name_to_id = {v: i for i, v in enumerate(vertex_names)}
        #: Pointer to the containing sets
        self._id = list(range(len(vertex_names)))
        #: Size of the set (_sz[_id[v]] is the size of the set that contains v)
        self._sz = [1] * len(vertex_names)

    def find(self, v):
        assert type(v) is int
        j = v

        while j != self._id[j]:
            self._id[j] = self._id[self._id[j]]
            j = self._id[j]

        return j

    def find_by_name(self, v_name):
        return self.find(self._name_to_id[v_name])

    def union_by_name(self, v_name, w_name):
        self.union(self.find_by_name(v_name), self.find_by_name(w_name))

    def union(self, v, w):
        assert type(v) is int
        assert type(w) is int
        i = self.find(v)
        j = self.find(w)

        if i == j:
            return

        if self._sz[i] < self._sz[j]:
            self._id[i] = j
            self._sz[j] += self._sz[i]

        else:
            self._id[j] = i

        self._sz[i] += self._sz[j]


@attr.s(frozen=True, auto_attribs=True)
class CopyNumberVariant:
    """Represent on CNV from one sample."""

    #: The chromosome.
    chrom: str
    #: The 0-based start position.
    pos_begin: int
    #: The 0-based end position.
    pos_end: int
    #: The kind of the CNV (del/dup).
    kind: str
    #: The source/caller of the CNV.
    source: str
    #: The sample this CNV was seen in.
    sample: str
    #: Annotation of the CNV.
    anno: typing.Dict[str, typing.Any]

    def recip_ovl(self, other: typing.TypeVar("CopyNumberVariant")) -> float:
        """Compute reciprocal overlap of self with other."""
        if self.chrom != other.chrom:
            return False
        elif self.pos_begin <= other.pos_end and other.pos_begin <= self.pos_end:
            a = min(self.pos_end, other.pos_end) - max(self.pos_begin, other.pos_begin)
            b = max(self.pos_end - self.pos_begin, other.pos_end - other.pos_begin)
            return a / b
        else:
            return 0.0


@attr.s(frozen=True, auto_attribs=True)
class ContigCnvs:
    """Store the CNVs for one contig with lookup table."""

    #: The contig name.
    contig: str
    #: The CopyNumberVariant objects.
    cnvs: typing.Tuple[CopyNumberVariant, ...]
    #: The interval lookup table.
    ncls: ncls.NCLS

    @staticmethod
    def from_cnvs(contig: str, cnvs: typing.Iterable[CopyNumberVariant]) -> typing.TypeVar(
        "ContigCnvs"
    ):
        """Build from name and list of CopyNumberVariant."""
        start = pd.Series([cnv.pos_begin for cnv in cnvs])
        ends = pd.Series([cnv.pos_end for cnv in cnvs])
        ids = pd.Series(range(0, len(cnvs)))
        lookup = ncls.NCLS(start.values, ends.values, ids.values)
        return ContigCnvs(contig, tuple(cnvs), lookup)


@attr.s(frozen=True, auto_attribs=True)
class CnvCluster:
    """Represent one cluster of CNVs."""

    #: The CopyNumberVariant objects.
    cnvs: typing.Tuple[CopyNumberVariant, ...]


def merge_headers(headers: typing.Iterable[vcfpy.Header]) -> vcfpy.Header:
    """Merge Headers for output."""
    res = None
    for header in headers:
        if res is None:
            res = vcfpy.Header(list(header.lines), samples=vcfpy.SamplesInfos(header.samples.names))
        else:
            for line in header.lines:
                if "ID" in getattr(line, "mapping", {}) and not res.has_header_line(
                    line.key, line.mapping["ID"]
                ):
                    res.add_line(line)
            res.samples.names += header.samples.names
    return res


def augment_header(header: vcfpy.Header) -> vcfpy.Header:
    if not header.has_header_line("INFO", "SVMETHOD"):
        header.add_info_line(
            {
                "ID": "SVMETHOD",
                "Number": 1,
                "Type": "String",
                "Description": "Type of approach used to detect SV",
            }
        )
    if not header.has_header_line("INFO", "SVTYPE"):
        header.add_info_line(
            {
                "ID": "SVTYPE",
                "Number": 1,
                "Type": "String",
                "Description": "Type of structural variant",
            }
        )
    if not header.has_header_line("INFO", "SVLEN"):
        header.add_info_line(
            {
                "ID": "SVLEN",
                "Number": ".",
                "Type": "Integer",
                "Description": "Difference in length between REF and ALT alleles",
            }
        )
    if not header.has_header_line("INFO", "CIPOS"):
        header.add_info_line(
            {
                "ID": "CIPOS",
                "Number": 2,
                "Type": "Integer",
                "Description": "Confidence interval around POS for imprecise variants",
            }
        )
    if not header.has_header_line("INFO", "CIEND"):
        header.add_info_line(
            {
                "ID": "CIEND",
                "Number": 2,
                "Type": "Integer",
                "Description": "Confidence interval around END for imprecise variants",
            }
        )
    return header


def process_contig_inner(contig: str, reader: vcfpy.Reader) -> typing.List[CopyNumberVariant]:
    try:
        _contig_iter = reader.fetch(contig)  # noqa: F841
    except ValueError as _e:  # noqa: F841
        return  # contig not in file, skip

    for record in reader.fetch(contig):
        if not record.ALT:
            logger.debug("Skipping %s (no CNV)", ";".join(record.ID))
        else:
            for sample, call in record.call_for_sample.items():
                yield CopyNumberVariant(
                    chrom=record.CHROM,
                    pos_begin=record.affected_start,
                    pos_end=record.INFO["END"],
                    kind=record.ALT[0].value,
                    source=SOURCE_GATK_GCNV,
                    sample=sample,
                    anno=call.data,
                )


def cluster_cnvs(contig_cnvs: ContigCnvs, min_ovl: float) -> typing.Tuple[CnvCluster]:
    logger.debug("overlapping for contig %s", contig_cnvs.contig)
    num_cnvs = len(contig_cnvs.cnvs)
    uf = UnionFind(range(num_cnvs))
    for i, cnv in enumerate(contig_cnvs.cnvs):
        for _start, _end, j in contig_cnvs.ncls.find_overlap(cnv.pos_begin, cnv.pos_end):
            other = contig_cnvs.cnvs[j]
            logger.debug("TEST: %s / %s / %f", cnv, other, cnv.recip_ovl(other))
            if i != j and cnv.kind == other.kind and cnv.recip_ovl(other) >= min_ovl:
                logger.debug("OVERLAP %s / %s", cnv, other)
                uf.union(i, j)

    out_ids = list(sorted(set(map(uf.find, range(num_cnvs)))))
    id_to_out = {v: k for k, v in enumerate(out_ids)}
    res = [list() for i in range(len(out_ids))]
    for i in range(num_cnvs):
        res[id_to_out[uf.find(i)]].append(contig_cnvs.cnvs[i])

    for value in res:
        logger.debug("RESULT (%d) / %s: %s", len(value), "|".join([v.sample for v in value]), value)

    return tuple([CnvCluster(tuple(x)) for x in res])


def cluster_to_record(cluster: CnvCluster, header: vcfpy.Header) -> vcfpy.Record:
    first = cluster.cnvs[0]
    # TODO: augment with start/end CI
    min_pos_begin = min((cnv.pos_begin + 1 for cnv in cluster.cnvs))
    mean_pos_begin = int(mean((cnv.pos_begin + 1 for cnv in cluster.cnvs)))
    max_pos_begin = max((cnv.pos_begin + 1 for cnv in cluster.cnvs))
    min_pos_end = min((cnv.pos_end for cnv in cluster.cnvs))
    mean_pos_end = int(min((cnv.pos_end for cnv in cluster.cnvs)))
    max_pos_end = max((cnv.pos_end for cnv in cluster.cnvs))

    info = {}
    if mean_pos_begin > max_pos_end:
        pos = (mean_pos_end + mean_pos_begin) / 2
        mean_pos_end = pos
        mean_pos_begin = pos

    info["END"] = mean_pos_end
    info["SVTYPE"] = first.kind
    if min_pos_end != max_pos_end:
        info["CIPOS"] = [min_pos_begin - mean_pos_begin, max_pos_begin - mean_pos_begin]
        info["CIEND"] = [min_pos_end - mean_pos_end, max_pos_end - mean_pos_end]
    if first.kind == "DEL":
        info["SVLEN"] = [-(mean_pos_end - mean_pos_begin)]
    else:  # first.kind == "DUP"
        info["SVLEN"] = [mean_pos_end - mean_pos_begin]
    sample_to_data = {}
    fmt = ["GT"]
    for cnv in cluster.cnvs:
        fmt = list(cnv.anno.keys())
        sample_to_data[cnv.sample] = cnv.anno
    calls = []
    for sample in header.samples.names:
        if sample in sample_to_data:
            calls.append(vcfpy.Call(sample, sample_to_data.get(sample)))
        else:
            calls.append(vcfpy.Call(sample, {"GT": "0/0", "CN": 2}))
    return vcfpy.Record(
        CHROM=first.chrom,
        POS=mean_pos_begin,
        ID="",
        REF="N",
        ALT=[vcfpy.SymbolicAllele(first.kind)],
        QUAL=None,
        FILTER=[],
        INFO=info,
        FORMAT=fmt,
        calls=calls,
    )


def process_contig(
    contig: str, readers: typing.Iterable[vcfpy.Reader], min_ovl: float, out_header: vcfpy.Header
):
    cnvs = []
    for reader in readers:
        source = reader.header._indices["source"][0].value
        if source not in SOURCE_MAP:
            raise Exception("Unknown source: %s" % source)
        source = SOURCE_MAP[source]
        logger.debug("File %s from source %s", os.path.basename(reader.path), source)
        cnvs += list(process_contig_inner(contig, reader))
    if not cnvs:
        return []  # no CNVs for contig
    else:
        logger.info(
            "Parsed a total of %d CNVs from %d VCFs for contig %s", len(cnvs), len(readers), contig
        )

    logger.info("Building contig CNV lookup...")
    contig_cnvs = ContigCnvs.from_cnvs(contig, cnvs)

    logger.info("Clustering CNVs...")
    cnv_clusters = cluster_cnvs(contig_cnvs, min_ovl)
    logger.info("Created %d clusters", len(cnv_clusters))

    logger.info("Converting clusters to records...")
    records = (cluster_to_record(cluster, out_header) for cluster in cnv_clusters)

    return sorted(
        records, key=lambda record: (record.POS, -record.INFO["END"], record.INFO["SVTYPE"])
    )


def run(args):
    logger.info("Starting gCNV VCF merging")
    logger.info("config = %s", args)

    with contextlib.ExitStack() as stack:
        logger.info("Open input files and merge headers...")
        readers = [stack.enter_context(vcfpy.Reader.from_path(path)) for path in args.in_vcf]
        out_header = augment_header(merge_headers([reader.header for reader in readers]))
        writer = stack.enter_context(vcfpy.Writer.from_path(args.out_vcf, out_header))

        logger.info("Processing contigs...")
        for contig_line in out_header.get_lines("contig"):
            records = process_contig(contig_line.mapping["ID"], readers, args.min_ovl, out_header)
            for record in records:
                writer.write_record(record)
        logger.info("Done processing contig %s.", contig_line.mapping["ID"])

    logger.info("All done. Have a nice day!")


def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument("out_vcf", metavar="OUT.vcf", help="Path to output VCF file")
    parser.add_argument("in_vcf", metavar="IN.vcf", nargs="+", help="Path to input VCF files ")
    parser.add_argument("--min-ovl", type=float, default=0.75, help="Minimal reciprocal overlap")
    parser.add_argument(
        "--verbose", "-v", default=False, action="store_true", help="Enable verbose mode"
    )

    args = parser.parse_args(argv)
    if args.verbose:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
