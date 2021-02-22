#!/usr/bin/env python3
"""Simple script that only keeps the first seen VCF header

Reads from STDIN, writes to STDOUT

Usage::

    $ cat INPUT.vcf | snappy-vcf_first_header >OUTPUT.vcf
"""

import sys

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def main():
    """Main entry point for the script"""
    seen_header = False
    for line in sys.stdin:
        is_header = line.startswith("#")
        is_chrom = line.startswith("#CHROM")
        if seen_header:
            if not is_header:
                print(line, file=sys.stdout, end="")
        else:  # not seen_header
            print(line, file=sys.stdout, end="")
        if is_chrom:
            seen_header = True


if __name__ == "__main__":
    sys.exit(main())
