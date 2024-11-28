#!/usr/bin/env python3
"""Generate genome regions from FAI file

Usage::

    $ snappy-genome_windows --fai-file FILE.fa.fai

    $ snappy-genome_windows --fai-file FILE.fa.fai --format bed --output-file OUT.bed
"""

import argparse
import csv
import fnmatch
import os
import re
import sys

from pathlib import Path


# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.genome_regions import GenomeRegion  # noqa: E402

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default for ``--window-size``
DEFAULT_WINDOW_SIZE = 10000000
#: Default for `--format``
DEFAULT_FORMAT = "regions"

#: Allowed values for ``--format``
CHOICES_FORMAT = ("regions", "bed")

#: Regular expression patterns to parse *.fai, *.genome, *.dict & fastq files
PATTERN_FAI = re.compile(r"^([^\s]+)\t([0-9]+)\t([0-9]+)\t([0-9]+)\t([0-9]+)\s*$")
PATTERN_GENOME = re.compile(r"^([^\s]+)\t([0-9]+)\s*$")
PATTERN_DICT = re.compile(r"^@SQ\tSN:([^\s]+)\tLN:([0-9]+).*$")
PATTERN_FASTA = re.compile(r"^\s*>\s*([^\s]+).*$")


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


def ignore_chroms(path_ref: str, ignored: set[str] = [], return_ignored: bool = False):
    path_ref = Path(path_ref).resolve()
    if Path(str(path_ref) + ".fai").exists():
        contigs = _parse_index(Path(str(path_ref) + ".fai"), PATTERN_FAI)
    elif Path(str(path_ref) + ".genome").exists():
        contigs = _parse_index(Path(str(path_ref) + ".genome"), PATTERN_GENOME)
    elif path_ref.with_suffix("dict").exists():
        contigs = _parse_index(path_ref.with_suffix("dict"), PATTERN_DICT, True)
    else:
        contigs = _read_fasta(path_ref)
    for contig_name, contig_length in contigs:
        m = matches_any(contig_name, ignored)
        if (m and return_ignored) or (not m and not return_ignored):
            yield contig_name, contig_length


def _parse_index(filename: Path, pattern: re.Pattern, allow_mismatch: bool = False):
    with open(filename, "rt") as f:
        for line in f:
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            m = pattern.match(line)
            if m:
                groups = m.groups()
                yield groups[0], int(groups[1])
            else:
                if not allow_mismatch:
                    raise ValueError(f"Unexpected record '{line}' in reference file '{filename}'")


def _read_fasta(filename: Path):
    contig_name = None
    contig_length = None
    with open(filename, "rt") as f:
        for line in f:
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            m = PATTERN_FASTA.match(line)
            if m:
                if contig_name:
                    yield contig_name, contig_length
                groups = m.groups()
                contig_name = groups[0]
                contig_length = 0
            else:
                contig_length += len(line)
    assert contig_name is not None, f"No contig found in reference file {filename}"
    yield contig_name, contig_length


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
