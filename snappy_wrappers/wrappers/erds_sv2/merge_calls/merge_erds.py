#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Merge ERDS call output files.
"""

import argparse
from collections import OrderedDict
import re
import sys
import typing

import vcfpy

from intervaltree import IntervalTree


class CnvCall(typing.NamedTuple):
    """Represent a CNV call"""

    #: SV type ('DEL'/'DUP')
    svtype: str
    #: Chromosome name.
    chrom: str
    #: 0-based start position
    begin: int
    #: 0-based end position
    end: int
    #: 0-based confidence interval around start position
    ci: typing.Tuple[int, int]
    #: 0-based confidence interval around end position
    ciend: typing.Tuple[int, int]

    def overlap(self, other):
        """Compute overlap between ``self`` and ``other``."""
        if self.chrom != other.chrom:
            return 0
        else:
            ovl_end = min(self.end, other.end)
            ovl_start = max(self.begin, other.begin)
            return max(ovl_end - ovl_start, 0)

    def length(self):
        """Return interval length."""
        return self.end - self.begin

    def reciprocal_overlap(self, other):
        """Compute reciprocal overlap for ``self``, and ``other``."""
        ovl = self.overlap(other)
        return min(ovl / self.length(), ovl / other.length())

    def merge(self, other):
        """Compute merge result of ``self`` and ``other``."""
        assert self.overlap(other) > 0
        assert self.svtype == other.svtype
        return CnvCall(
            svtype=self.svtype,
            chrom=self.chrom,
            begin=self.begin + (other.begin - other.begin) // 2,
            end=self.end + (other.end - other.end) // 2,
            ci=(min(self.ci[0], other.ci[0]), max(self.ci[1], other.ci[1])),
            ciend=(min(self.ciend[0], other.ciend[0]), max(self.ciend[1], other.ciend[1])),
        )


def process_input(path_in, region, ovl_thresh, tree, call_list):
    with vcfpy.Reader.from_path(path_in) as reader:
        print("Fetching region {} for {}".format(region, path_in), file=sys.stderr)
        fix_reader(reader)
        chrom, itv = region.split(":")
        begin, end = list(map(int, map(lambda x: x.replace(",", ""), itv.split("-")[:2])))
        try:
            iter = reader.fetch(chrom, begin - 1, end)
        except ValueError as e:
            print("Problem fetching region: {}".format(e), file=sys.stderr)
            iter = []
        # Collect CNV calls
        calls = []
        for record in iter:
            # print(record.CHROM, record.POS)
            if record.INFO.get("SVTYPE") not in ("DEL", "DUP"):
                print(
                    "Skipping record with SVTYPE {}".format(record.INFO.get("SVTYPE")),
                    file=sys.stderr,
                )
                continue
            calls.append(
                CnvCall(
                    svtype=record.INFO["SVTYPE"],
                    chrom=record.CHROM,
                    begin=record.affected_start,
                    end=record.INFO.get("END"),
                    ci=(-1, 1),
                    ciend=(-1, 1),
                )
            )
        # Merge with previous calls (with all overlapping):
        for call in calls:
            matches = tree.search(call.begin, call.end)
            if not matches:
                idx = len(call_list)
                call_list.append(call)
                tree.addi(call.begin + call.ci[0], call.end + call.ciend[1], idx)
            else:
                good = []
                for itv in matches:
                    other = call_list[itv.data]
                    ovl = other.reciprocal_overlap(call)
                    if ovl >= ovl_thresh and other.svtype == call.svtype:
                        good.append(itv)
                for g in good:
                    tree.remove(g)
                    merged = call_list[g.data].merge(call)
                    call_list[g.data] = merged
                    tree.addi(merged.begin + merged.ci[0], merged.end + merged.ciend[1], g.data)


def write_region(tree, call_list, writer):
    calls = [call_list[itv.data] for itv in tree.items()]
    calls = sorted(calls, key=lambda x: x.begin)
    for call in calls:
        if call.svtype == "DEL":
            svlen = -call.length()
        else:
            svlen = call.length()
        infos = OrderedDict(
            [
                ("SVTYPE", call.svtype),
                ("END", call.end),
                ("CI", call.ci),
                ("CIEND", call.ciend),
                ("SVLEN", [svlen]),
            ]
        )
        record = vcfpy.Record(
            call.chrom,
            call.begin + 1,
            [],
            "N",
            [vcfpy.SymbolicAllele(call.svtype)],
            None,
            [],
            infos,
            [],
            OrderedDict(),
        )
        writer.write_record(record)


def full_chromosomes(fai_path):
    """Return list of regions of all chromosomes of VCF reader."""
    with open(fai_path, "rt") as fai:
        for line in fai:
            line = line.strip()
            if line:
                arr = line.split("\t")
                chrom, length = arr[:2]
                yield "{}:1-{}".format(chrom, length)


def fix_reader(reader):
    """Fix header of reader."""
    reader.header.add_filter_line(
        OrderedDict([("ID", "PASS"), ("Description", "All filters passed")])
    )
    reader.header.add_info_line(
        OrderedDict(
            [
                ("ID", "PRECISE"),
                ("Number", 0),
                ("Type", "Flag"),
                ("Description", "Precise structural variation"),
            ]
        )
    )


def run(args):
    # Use header from first ERDS output file.
    with vcfpy.Reader.from_path(args.inputs[0]) as reader:
        fix_reader(reader)
        header = vcfpy.header_without_lines(reader.header.copy(), (("contig", None)))
        header.samples.names = []
    header.add_info_line(
        OrderedDict(
            [
                ("ID", "CI"),
                ("Number", 2),
                ("Type", "Integer"),
                ("Description", "Confidence interval around POS"),
            ]
        )
    )
    header.add_info_line(
        OrderedDict(
            [
                ("ID", "CIEND"),
                ("Number", 2),
                ("Type", "Integer"),
                ("Description", "Confidence interval around END"),
            ]
        )
    )
    # Add contigs from VCF file.
    with open(args.fai, "rt") as fai:
        for line in fai:
            line = line.strip()
            if line:
                contig, length = line.split("\t")[:2]
                header.add_contig_line(OrderedDict([("ID", contig), ("length", length)]))
    # Add command line header
    header.add_line(vcfpy.HeaderLine(key="merge_erds_cmd", value=" ".join(sys.argv)))
    # Get list of regions to process (from args or VCF header).
    if args.regions:
        regions = args.regions
    else:
        regions = []
        for region in full_chromosomes(args.fai):
            contig = region.split(":")[0]
            if re.match(args.chrom_regex, contig):
                regions.append(region)
    # Create writer and iterate over all regions from all input files.
    with vcfpy.Writer.from_path(args.output, header) as writer:
        for region in regions:
            tree = IntervalTree()
            call_list = []
            for path_in in args.inputs:
                process_input(path_in, region, args.overlap, tree, call_list)
            write_region(tree, call_list, writer)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Merging of ERDS output files")

    group = parser.add_argument_group("General Options")
    group.add_argument("-v", "--verbose", default=0, action="count")

    group = parser.add_argument_group("Input / Output Options")
    group.add_argument("--chrom-regex", default=r"^(chr)?([0-9XY]+|MT?)$")
    group.add_argument(
        "--input",
        nargs="+",
        action="append",
        default=[],
        dest="inputs",
        required=True,
        help="input VCF file",
    )
    group.add_argument("--fai", required=True, help="Path to FAI file of reference")
    group.add_argument("--output", help="output VCF file", required=True)
    group.add_argument(
        "--region",
        type=str,
        required=False,
        default=[],
        action="append",
        dest="regions",
        nargs="+",
        help=("region(s) to limit analysis to"),
    )

    group = parser.add_argument_group("Merging Options")
    group.add_argument("--overlap", default=0.8, type=str, help="Required reciprocal overlap")

    args = parser.parse_args(argv)
    args.inputs = [i for lst in args.inputs for i in lst]  # flatten
    args.regions = [r for lst in args.regions for r in lst]  # flatten

    print("Arguments: {}".format(args), file=sys.stderr)

    return run(args)


if __name__ == "__main__":
    sys.exit(main())
