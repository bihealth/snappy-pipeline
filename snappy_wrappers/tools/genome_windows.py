#!/usr/bin/env python3
"""The contigs of a FAI file, without the ignored ones"""

import csv
import re

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def matches_any(query, patterns):
    """Whether any regular expression matches ``query`` (``re.search``, i.e. unanchored)."""
    for pattern in patterns:
        if re.search(pattern, query):
            return True
    return False


def yield_contigs(fai_file, ignore_contigs=None):
    """Yield contig names."""
    csv_reader = csv.reader(fai_file, delimiter="\t")
    for record in csv_reader:
        chrom, _, _, _, _ = record
        if not matches_any(chrom, ignore_contigs or []):
            yield chrom
