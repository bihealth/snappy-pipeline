#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datetime
import logging

import sys

import vcfpy


def run(args):
    logging.info("Moving FILTER to INFO/ORIG_FILTER")
    logger = logging.getLogger("")  # root logger
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    time_start = datetime.datetime.now()

    with vcfpy.Reader.from_path(args.input_vcf) as reader:
        header = reader.header
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "ORIG_FILTER"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "Original FILTER value"),
                ]
            )
        )
        with vcfpy.Writer.from_path(args.output_vcf, header) as writer:
            for record in reader:
                if record.FILTER and record.FILTER != ["PASS"]:
                    record.INFO["ORIG_FILTER"] = record.FILTER
                    record.FILTER = []
                writer.write_record(record)

    time_end = datetime.datetime.now()
    spent = time_end - time_start
    logging.info("Spent %.1f s", spent.total_seconds())


def main(argv=None):
    """Program's main entry point (before parsing command line arguments)."""

    # Setup command line parser
    parser = argparse.ArgumentParser(description="Move FILTER column to INFO")

    parser.add_argument(
        "--verbose",
        "-v",
        dest="verbose",
        default=False,
        action="store_true",
        help="Enable verbose logging",
    )

    group = parser.add_argument_group("Input / Output")
    group.add_argument("--input-vcf", required=True, help="Path to input VCF file")
    group.add_argument("--output-vcf", required=True, help="Path to output VCF file")

    # Parse arguments, postprocess and kick-off program.
    args = parser.parse_args(argv)

    run(args)


if __name__ == "__main__":
    sys.exit(main())
