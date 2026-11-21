# -*- coding: utf-8 -*-
"""Code for genome regions"""


class GenomeRegion:
    """Genome region with half-open intervals"""

    def __init__(self, chrom, begin, end):
        self.chrom = chrom
        self.begin = begin
        self.end = end

    def as_bed(self):
        """Return BED reprentation (half-open intervals)"""
        tpl = "{}\t{}\t{}"
        return tpl.format(self.chrom, self.begin, self.end)

    def human_readable(self, with_comma=True):
        """Return human readable string"""
        if with_comma:
            tpl = "{}:{:,}-{:,}"
        else:
            tpl = "{}:{:}-{:}"
        return tpl.format(self.chrom, self.begin + 1, self.end)
