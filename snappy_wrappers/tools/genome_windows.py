#!/usr/bin/env python3
"""Generate genome regions from FAI file

Usage::

    $ snappy-genome_windows --fai-file FILE.fa.fai

    $ snappy-genome_windows --fai-file FILE.fa.fai --format bed --output-file OUT.bed
"""

import argparse
import csv
import fnmatch
import sys

from snappy_wrappers.genome_regions import GenomeRegion

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default for ``--window-size``
DEFAULT_WINDOW_SIZE = 10000000
#: Default for `--format``
DEFAULT_FORMAT = "regions"

#: Allowed values for ``--format``
CHOICES_FORMAT = ("regions", "bed")


def matches_any(query, patterns):
    for pattern in patterns:
        if fnmatch.fnmatch(query, pattern):
            return True
    return False


def yield_contigs(fai_file, ignore_chroms=None):
    """Yield contig names."""
    csv_reader = csv.reader(fai_file, delimiter="\t")
    for record in csv_reader:
        chrom, _, _, _, _ = record
        if not matches_any(chrom, ignore_chroms or []):
            yield chrom


def yield_regions(fai_file, window_size, subtract_end=0, ignore_chroms=None, padding=0):
    """Yield GenomeRegion by GenomeRegion

    ``subtract_end`` -- for 0/1 based coordinates.
    ``padding`` -- padding of windows towards ech size
    """
    csv_reader = csv.reader(fai_file, delimiter="\t")
    for record in csv_reader:
        chrom, length, _, _, _ = record
        if matches_any(chrom, ignore_chroms or []):
            continue  # skip this chromosomes
        length = int(length)
        # process chromosome
        begin = 0
        while begin < length:
            end = begin + window_size
            if end > length:
                end = length
            # Compute padding without reaching over contig
            pad_left = min(padding, begin)
            pad_right = min(padding, length - end)
            # Construct GenomeRegion object
            region = GenomeRegion(chrom, begin - pad_left, end + pad_right)
            region.end -= subtract_end
            if begin < end - subtract_end:
                yield region
            begin = end


def run(args):
    """Main entry point after parsing command line arguments"""
    yielded = 0
    for region in yield_regions(
        args.fai_file, args.window_size, args.subtract_end, args.ignore_chroms
    ):
        if args.format == "regions":
            print(region.human_readable(), file=args.output_file)
        else:  # args.format == 'bed'
            print(region.as_bed(), file=args.output_file)
        yielded += 1
        if args.count and yielded >= args.count:
            break


def create_parser():
    """Construct and return the command line parser"""
    parser = argparse.ArgumentParser(description="Build genome windows from FAI file")
    parser.add_argument(
        "--fai-file",
        type=argparse.FileType("rt"),
        required=True,
        help="FAI file to generate windows for",
    )
    parser.add_argument(
        "--window-size", type=int, default=DEFAULT_WINDOW_SIZE, help="Window length to generate"
    )
    parser.add_argument(
        "--output-file",
        type=argparse.FileType("wt"),
        default=sys.stdout,
        help="Output file for regions",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=CHOICES_FORMAT,
        default=DEFAULT_FORMAT,
        help="Use region strings (defaults) or BED format",
    )
    parser.add_argument(
        "--subtract-end",
        type=int,
        default=0,
        help="Optionally, strip one from end, sometimes required for BED",
    )
    parser.add_argument(
        "--ignore-chroms",
        default=[],
        action="append",
        nargs="+",
        help="Patterns for contigs to ignore",
    )
    parser.add_argument(
        "--count", default=0, type=int, help="Number of windows to limit to, if any"
    )
    return parser


def main(argv=None):
    """Main entry point, includes parsing of command line arguments"""
    parser = create_parser()
    args = parser.parse_args(argv)
    args.ignore_chroms = [item for sublist in args.ignore_chroms for item in sublist]
    run(args)


if __name__ == "__main__":
    sys.exit(main())
