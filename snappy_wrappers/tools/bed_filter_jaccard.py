#!/usr/bin/env python3
"""Jaccard-index based filter for BED files

The default configuration works for two BED3 files.  Otherwise, you have to adjust the
``--num-cols-first`` and ``--num-cols-second`` arguments.

Usage::

    bedtools intersect -wao -a FIRST.bed -b SECOND.bed \
    | filter_jaccard.py --threshold 0.8 \
    > out.bed
"""

import argparse
import csv
import sys

#: Default value for ``--operation``
DEFAULT_OPERATION = "subtract"
#: Default value for ``--threshold``
DEFAULT_THRESHOLD = 0.8
#: Default value for ``--num-cols-first``
DEFAULT_NUM_COLS_FIRST = 3
#: Default value for ``--num-cols-second``
DEFAULT_NUM_COLS_SECOND = 3

#: Choices for ``--operation``
CHOICES_OPERATION = ("subtract", "intersect")


class Stats:
    """Simple counter for above/not about threshold"""

    def __init__(self, row, yes=0, no=0):
        self.row = row
        self.yes = yes
        self.no = no


def run(args):
    """Program entry point after parsing the command line"""
    stats = {}
    # Read in bedtools intersect results and compute Jaccard indices
    reader = csv.reader(args.input_file, delimiter="\t")
    for row in reader:
        chrom1, begin1, end1 = row[:3]
        _, begin2, end2 = row[args.num_cols_first : args.num_cols_first + 3]
        intersection = row[args.num_cols_first + args.num_cols_second]
        begin1, end1, begin2, end2, intersection = (
            int(begin1),
            int(end1),
            int(begin2),
            int(end2),
            int(intersection),
        )
        union = max(end1, end2) - min(begin1, begin2)
        jaccard = intersection / union
        stats.setdefault((chrom1, begin1, end1), Stats(row))
        if jaccard >= args.threshold:
            stats[(chrom1, begin1, end1)].yes += 1
        else:
            stats[(chrom1, begin1, end1)].no += 1
    # Write out resulting intervals
    for _, stats in sorted(stats.items()):
        if (args.operation == "subtract" and not stats.yes) or (
            args.operation == "intersect" and stats.yes
        ):
            print("\t".join(stats.row), file=args.output_file)


def main(argv=None):
    """Program entry point for parsing command line"""
    parser = argparse.ArgumentParser(
        description="Compute BED file overlaps using the Jaccard index"
    )
    parser.add_argument(
        "--input-file",
        type=argparse.FileType("rt"),
        default=sys.stdin,
        help="bedtools intersect file",
    )
    parser.add_argument(
        "--output-file", type=argparse.FileType("wt"), default=sys.stdout, help="Resulting BED file"
    )
    parser.add_argument(
        "--operation",
        type=str,
        default=DEFAULT_OPERATION,
        choices=CHOICES_OPERATION,
        help="Operation to perform; default: %s" % DEFAULT_OPERATION,
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Treshold on required Jaccard index; default: %f" % DEFAULT_THRESHOLD,
    )
    parser.add_argument(
        "--num-cols-first",
        type=int,
        default=DEFAULT_NUM_COLS_FIRST,
        help=(
            "Number of columns in first original file (to be filtered down); default: %d"
            % DEFAULT_NUM_COLS_FIRST
        ),
    )
    parser.add_argument(
        "--num-cols-second",
        type=int,
        default=DEFAULT_NUM_COLS_SECOND,
        help=(
            "Number of columns in second original file (to be used for filtering down); "
            "default: %d" % DEFAULT_NUM_COLS_SECOND
        ),
    )

    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    sys.exit(main())
