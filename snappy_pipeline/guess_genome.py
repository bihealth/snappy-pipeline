# -*- coding: utf-8 -*-
"""Helper code for guessing UCSC genome IDs from reference paths
"""

from collections import OrderedDict

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Tokens to look for in the path to the reference for guessing UCSC genome
#: ID.  This is a projection to the UCSC ids only.
TOKEN_TO_UCSC_ID = OrderedDict(
    [
        ("GRCh37", "hg19"),
        ("hg19", "hg19"),
        ("GRCh38", "hg38"),
        ("hg38", "hg38"),
        ("NCBIM37", "mm9"),
        ("mm9", "mm9"),
    ]
)


def guess_ucsc_genome_id(reference):
    """Try to guess UCSC genome ID from file system path ``reference``

    The UCSC genome ID could be ``hg19``, ``hg38`` etc.
    """
    for token, result in TOKEN_TO_UCSC_ID.items():
        if token in reference:
            return result
    return "hg19"  # default
