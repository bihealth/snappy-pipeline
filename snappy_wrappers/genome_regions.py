# -*- coding: utf-8 -*-
"""Code for genome regions"""

import re


class GenomeRegion:
    """Genome region with half-open intervals"""

    @staticmethod
    def from_bed_line(line):
        """Return GenomeRegion from BED line"""
        chrom, begin, end = line.strip().split("\t")
        return GenomeRegion(chrom, int(begin), int(end))

    @staticmethod
    def from_human_readable(human_readable):
        """Parse human-readable genome description into ``GenomeRegion``"""
        human_readable = human_readable.replace(",", "")
        m = re.match("^(?P<chrom>.*?):(?P<start>[0-9]+)-(?P<end>[0-9]+)$", human_readable)
        if not m:
            raise ValueError("Invalid region string: {}".format(human_readable))
        else:
            return GenomeRegion(m.group("chrom"), int(m.group("start")) - 1, int(m.group("end")))

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

    @property
    def length(self):
        """Return length"""
        return self.end - self.begin

    def overlaps(self, other):
        """Return whether the region overlas with ``other``"""
        if self.chrom != other.chrom:
            return False
        if self.begin < other.end and other.begin < self.end:
            return True
        else:
            return False

    def extend(self, by):
        """Extend genome region at both sides by ``by``"""
        return GenomeRegion(self.chrom, max(0, self.begin - by), self.end + by)

    def jaccard(self, other):
        """Return jaccard index between ``self`` and ``other``"""
        if self.chrom != other.chrom:
            return 0
        pos_begin = max(self.begin, other.begin)
        pos_end = min(self.end, other.end)
        if pos_end >= pos_begin:
            overlap = pos_end - pos_begin
        else:
            overlap = 0
        union = self.length + other.length - overlap
        return overlap / union

    def __str__(self):
        tpl = "GenomeRegion({}, {}, {})"
        return tpl.format(*list(map(repr, [self.chrom, self.begin, self.end])))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if self is other:
            return True
        else:
            return (self.chrom, self.begin, self.end) == (other.chrom, other.begin, other.end)

    def __neq__(self, other):
        return not self.__eq__(other)
