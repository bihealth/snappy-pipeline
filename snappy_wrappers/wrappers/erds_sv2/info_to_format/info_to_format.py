#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Move the SV2 INFO fields that belong to FORMAT where they belong.
"""

import argparse
import itertools
import logging
import sys

import vcfpy

# White-listed chromosomes.
CHROMS = tuple(itertools.chain(map(str, range(1, 23)), ("X", "Y")))


def full_chromosomes(reader):
    """Return list of regions of all chromosomes of VCF reader."""
    for line in reader.header.get_lines("contig"):
        if line.id in CHROMS:
            name = line.id
            length = line.length
            yield "{}:{}-{}".format(name, 1, length)


class InfoToFormatApp:
    """Conversion from INFO to FORMAT field"""

    def __init__(self, args):
        #: Command line arguments.
        self.args = args
        # Setup the logging.
        self._setup_logging()

    def _setup_logging(self):
        logging.basicConfig(
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
        )
        logger = logging.getLogger("")
        if self.args.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)

    def _print_header(self):
        logging.info("INFO to FORMAT")
        logging.info("Arguments: %s", self.args)

    def _process_region(self, region, writer):
        """Process a single region and write its result to the writer."""

    def _augment_header(self, header):
        """Augment header information"""
        header = self._augment_filter(header)
        header = self._augment_info(header)
        header = self._augment_format(header)
        return header

    def _augment_filter(self, header):
        """Augment header for FILTER column"""
        return header

    def _augment_info(self, header):
        """Augment header for INFO column"""
        return header

    def _augment_format(self, header):
        """Augment header for FORMAT column"""
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SVMETHOD"),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", "Type of approach used to detect SV"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "FT"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "Filters from FILTER column of genotyping"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "DFT"),
                    ("Number", "1"),
                    ("Type", "String"),
                    (
                        "Description",
                        "Stringent filter status, recommended for de novo mutation discovery",
                    ),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "RG"),
                    ("Number", "1"),
                    ("Type", "Float"),
                    ("Description", "Median Phred-adjusted REF genotype likelihood"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "AF"),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", "Alternate allele frequency,in the range (0,1)"),
                ]
            )
        )
        return header

    def run(self):
        self._print_header()
        with vcfpy.Reader.from_path(self.args.input) as reader:
            # If no regions are given, fall back to all chromosomes.
            regions = full_chromosomes(reader)
            # Extend header with new lines.
            header = self._augment_header(reader.header)
            # Open the VCF writer for writing and process each region.
            with vcfpy.Writer.from_path(self.args.output, header) as writer:
                for region in regions:
                    logging.info("Processing %s", region)
                    try:
                        records = reader.fetch(region)
                    except ValueError:
                        records = []
                        logging.warning("Could not fetch records for %s", region)
                    for record in records:
                        record = self._process(record)
                        writer.write_record(record)

    def _process(self, record):
        """Process ``record``."""
        filters = list(record.FILTER)
        if "PASS" in filters and len(filters) > 0:
            filters = [f for f in filters if f != "PASS"]
        elif not filters:
            filters = ["PASS"]
        record.INFO["SVMETHOD"] = self.args.svmethod
        record.add_format("FT", filters)
        record.add_format("RG", record.INFO.get("REF_GTL"))
        record.add_format("DFT", record.INFO.get("DENOVO_FILTER"))
        record.add_format("AF", record.INFO.get("AF"))
        record.FILTER = []
        del record.INFO["REF_GTL"]
        del record.INFO["DENOVO_FILTER"]
        del record.INFO["AF"]
        return record


def main(argv=None):
    parser = argparse.ArgumentParser(description="Move SV2 INFO to FORMAT fields")

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("General Options")
    group.add_argument("-v", "--verbose", default=0, action="count")

    group = parser.add_argument_group("Input / Output Options")
    group.add_argument("--svmethod", required=True, help="name to put into SVMETHOD INFO field")
    group.add_argument("--input", required=True, help="input VCF file")
    group.add_argument("--output", help="output VCF file", default="/dev/stdout")

    args = parser.parse_args(argv)
    return InfoToFormatApp(args).run()


if __name__ == "__main__":
    sys.exit(main())
