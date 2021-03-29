#!/usr/bin/env python3
"""Helper tool to convert CNVetti coverage output to hom. DEL calls.
"""

import argparse
import contextlib
import datetime
import logging
import sys
import typing

import attr
import logzero
import vcfpy
from logzero import logger


@attr.s(frozen=True, auto_attribs=True)
class HomDel:
    """Represent of hom. DEL from one sample."""

    #: The chromosome.
    chrom: str
    #: The 0-based start position.
    pos_begin: int
    #: The 0-based end position.
    pos_end: int
    #: The raw read count
    raw_reads: int
    #: The total number of target bases
    target_bases: int
    #: The name of the sample
    sample: str

    def extend(self, other: typing.TypeVar("HomDel")) -> typing.TypeVar("HomDel"):
        assert self.chrom == other.chrom
        return HomDel(
            self.chrom,
            min(self.pos_begin, other.pos_begin),
            max(self.pos_end, other.pos_end),
            self.raw_reads + other.raw_reads,
            self.target_bases + other.target_bases,
            other.sample,
        )

    @staticmethod
    def from_record(record: vcfpy.Record) -> typing.TypeVar("HomDel"):
        call = record.calls[0]
        return HomDel(
            record.CHROM,
            record.POS,
            record.INFO["END"],
            call.data["RCV"],
            record.INFO["END"] - record.POS + 1,
            call.sample,
        )

    def to_record(self) -> vcfpy.Record:
        return vcfpy.Record(
            self.chrom,
            self.pos_begin,
            [],
            "N",
            [vcfpy.SymbolicAllele("DEL")],
            None,
            [],
            {
                "END": self.pos_end,
                "SVLEN": [self.pos_end - self.pos_begin + 1],
                "SVMETHOD": "cnvetti-homdel-0.2",
            },
            ["GT", "CN", "RCV", "LCV"],
            [
                vcfpy.Call(
                    self.sample,
                    {
                        "GT": "1",
                        "CN": 0.0,
                        "RCV": self.raw_reads,
                        "LCV": self.raw_reads / self.target_bases,
                    },
                )
            ],
        )


def build_header(header_in: vcfpy.Header) -> vcfpy.Header:
    result = vcfpy.Header(
        lines=[
            vcfpy.HeaderLine(key="fileformat", value="VCFv4.2"),
            vcfpy.HeaderLine(key="fileDate", value=datetime.datetime.now().strftime(r"%Y%m%d")),
            vcfpy.HeaderLine(key="source", value="CNVetti::homdel"),
            vcfpy.AltAlleleHeaderLine.from_mapping(
                {
                    "ID": "DEL",
                    "Description": "The record describes a deletion (decrease in coverage)",
                }
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "END",
                    "Number": 1,
                    "Type": "Integer",
                    "Description": "End position of the variant described in this record",
                }
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "SVTYPE",
                    "Number": 1,
                    "Type": "String",
                    "Description": "Type of structural variant",
                }
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "SVLEN",
                    "Number": ".",
                    "Type": "Integer",
                    "Description": "Difference in length between REF and ALT alleles",
                }
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                {
                    "ID": "SVMETHOD",
                    "Number": 1,
                    "Type": "String",
                    "Description": "Type of approach used to detect SV",
                }
            ),
            vcfpy.FormatHeaderLine.from_mapping(
                {"ID": "GT", "Number": 1, "Type": "String", "Description": "Genotype"}
            ),
            vcfpy.FormatHeaderLine.from_mapping(
                {
                    "ID": "CN",
                    "Number": 1,
                    "Type": "Float",
                    "Description": "Copy number of the copy number variant",
                }
            ),
            vcfpy.FormatHeaderLine.from_mapping(
                {"ID": "RCV", "Number": 1, "Type": "Float", "Description": "Raw coverage value"}
            ),
            vcfpy.FormatHeaderLine.from_mapping(
                {
                    "ID": "LCV",
                    "Number": 1,
                    "Type": "Float",
                    "Description": "Length-normalized coverage value",
                }
            ),
        ],
        samples=vcfpy.SamplesInfos(header_in.samples.names),
    )
    for line in header_in.lines:
        if line.key == "contig":
            result.add_contig_line({"ID": line.id, "length": line.length})
    return result


def process_contig(
    contig: str, reader: vcfpy.Reader, out_header: vcfpy.Header, max_rcv: float, max_lcv: float
):
    cnvs = []
    curr = None
    for record in reader:
        call = record.calls[0]
        if call.data["RCV"] < max_rcv or call.data["LCV"] < max_lcv:
            if curr:
                curr = curr.extend(HomDel.from_record(record))
            else:
                curr = HomDel.from_record(record)
        elif curr:
            yield curr.to_record()
            curr = None
    if curr:
        yield curr.to_record()


def run(args):
    logger.info("Starting to convert coverage to hom dels")
    logger.info("config = %s", args)

    with contextlib.ExitStack() as stack:
        logger.info("Open input and output file...")
        reader = stack.enter_context(vcfpy.Reader.from_path(args.in_vcf))
        out_header = build_header(reader.header)
        writer = stack.enter_context(vcfpy.Writer.from_path(args.out_vcf, out_header))

        logger.info("Processing contigs...")
        for contig_line in out_header.get_lines("contig"):
            for record in process_contig(
                contig_line.mapping["ID"], reader, out_header, args.max_rcv, args.max_lcv
            ):
                writer.write_record(record)
            logger.info("Done processing contig %s.", contig_line.mapping["ID"])

    logger.info("All done. Have a nice day!")


def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument("out_vcf", metavar="OUT.vcf", help="Path to output VCF file")
    parser.add_argument("in_vcf", metavar="IN.vcf", help="Path to input VCF files ")
    parser.add_argument(
        "--max-rcv", type=int, default=10, help="Maximum number of raw fragment count"
    )
    parser.add_argument(
        "--max-lcv", type=float, default=0.02, help="Maximum number of length-normalized fragments"
    )
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
