#!/usr/bin/env python3
"""Fix common issues with VCF files.

At the moment, these issues are:

- If no contig lines are present then they are added just before writing out
  the '^#CHROM...' line.
- If the file is empty then a minimal header is written out.
"""

import argparse
from functools import reduce
import gzip
import sys

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: No line has been read yet
STATE_INITIAL = 0
#: At least one header line has been read
STATE_SEEN_LINE = 1
#: Header has been read
STATE_HEADER_DONE = 2

PATTERN_CONTIG = "##contig="
PATTERN_CHROM = "#CHROM"


class VcfFileFixer:
    """Implementation of file fixing.

    Care has been taken to only look at lines when necessary as doing this in
    Python quickly becomes an I/O bottleneck.
    """

    def __init__(self, args, argv):
        #: Arguments as parsed from command line.
        self.args = args
        #: Command line arguments
        self.argv = tuple(argv or [])
        #: Output VCF file
        self.output_vcf = self.args.output_vcf
        #: Current state
        self.state = STATE_INITIAL
        #: Whether or not we have seen a contig line
        self.seen_contig_line = False
        #: State handlers
        self.handlers = (
            self._run_initial,  # 0
            self._run_seen_line,  # 1
            self._run_header_done,  # 2
        )
        #: Load contig lengths if path to FAIDX is given.
        self.contig_lengths = self._load_contig_lengths()

    def _load_contig_lengths(self):
        """Return contig lengths (tuple of pairs, so order is kept) or None"""
        if not self.args.faidx:
            return tuple()
        result = []
        for line in self.args.faidx:
            arr = line.strip().split("\t")
            result.append((arr[0], int(arr[1])))
        return tuple(result)

    def run(self):
        """Run the processing, simply dispatched to handlers."""
        for line in self.args.input_vcf:
            self.handlers[self.state](line)
        if self.state == STATE_INITIAL:
            self._print_dummy_header()

    def _run_initial(self, line):
        """Handle line processing in state `STATE_INITIAL```."""
        self.state = STATE_SEEN_LINE
        if line.startswith(PATTERN_CONTIG):
            self.seen_contig_line = True
        elif line.startswith(PATTERN_CHROM):
            self._handle_last_line()  # insert contig lines if necessary
        print(line, end="", file=self.output_vcf)

    def _run_seen_line(self, line):
        if line.startswith(PATTERN_CONTIG):
            self.seen_contig_line = True
        elif line.startswith(PATTERN_CHROM):
            self._handle_last_line()  # insert contig lines if necessary
        print(line, end="", file=self.output_vcf)

    def _run_header_done(self, line):
        """Header reading done, simply write to output line by line."""
        print(line, end="", file=self.output_vcf)

    def _handle_last_line(self):
        """Handler for transition to ``STATE_HEADER_DONE``."""
        self.state = STATE_HEADER_DONE
        if not self.seen_contig_line:
            self._print_contig_lines()
        self._print_fix_vcf_line()

    def _print_contig_lines(self):
        for name, length in self.contig_lengths:
            tpl = "##contig=<ID={id_},length={length}>"
            print(tpl.format(id_=name, length=length), file=self.output_vcf)

    def _print_dummy_header(self):
        """Print dummy header so the file is not empty."""
        print("##fileformat=VCFv4.3", file=self.output_vcf)
        self._print_contig_lines()
        self._print_fix_vcf_line()
        print("##fixVcfNote=The input file was empty, this is a dummy header")
        if self.args.samples:
            suff = "\tFORMAT\t{}".format("\t".join(self.args.samples))
        else:
            suff = ""
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{}".format(suff), file=self.output_vcf)

    def _print_fix_vcf_line(self):
        print("##fixVcf={}".format(" ".join(self.argv)))


def run(args, argv):
    """Run the filtration"""
    VcfFileFixer(args, argv).run()


def main(argv=None):
    """Program's main entry point"""
    parser = argparse.ArgumentParser(
        description="Apply various soft-filters on the WGS SV VCF file"
    )

    group = parser.add_argument_group("General Options")
    group.add_argument(
        "--input-vcf", type=argparse.FileType("rt"), required=True, help="input VCF file"
    )
    group.add_argument(
        "--output-vcf", type=argparse.FileType("wt"), help="output VCF file", default=sys.stdout
    )
    group.add_argument(
        "--faidx", type=argparse.FileType("rt"), help="FAI file for generating ##contig lines"
    ),
    group.add_argument("--sample", nargs="+", default=[], action="append", dest="samples")

    args = parser.parse_args(argv)
    args.samples = reduce(lambda x, y: x + y, args.samples)

    # also support compressed files (gzip or bgzip)
    if args.input_vcf.name.endswith(".gz"):
        args.input_vcf = gzip.open(args.input_vcf.name, "rt")

    run(args, argv or sys.argv)


if __name__ == "__main__":
    sys.exit(main())
