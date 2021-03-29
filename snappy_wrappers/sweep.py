"""Helpers for sweeping algorithms"""

import bisect

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class GenomePosition:
    """Represents a genome position with chromosome.

    Note that only positions from one chromosome are totally ordered.  Between chromosomes, there
    is no ordering.
    """

    LARGE_NUMBER = 10 ** 30

    def __init__(self, chrom, pos):
        self.chrom = chrom
        self.pos = pos

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return False
        return self.pos < other.pos

    def __gt__(self, other):
        if self.chrom != other.chrom:
            return False
        return other > self

    def __eq__(self, other):
        if self.chrom != other.chrom:
            return False
        return (self.chrom, self.pos) == (other.chrom, other.pos)

    def __le__(self, other):
        if self.chrom != other.chrom:
            return False
        return not other < self

    def __ge__(self, other):
        if self.chrom != other.chrom:
            return False
        return not self < other

    def delta(self, other):
        """Return offset between positions or LARGE_NUMBER if on different
        chromosomes
        """
        if self.chrom != other.chrom:
            return self.LARGE_NUMBER
        return self.pos - other.pos

    def in_range(self, chrom, begin, end):
        """Return true if within range on chromosome"""
        # TODO(holtgrewe): Use genomic range?
        if self.chrom != chrom:
            return False
        return self.pos >= begin and self.pos < end

    def __str__(self):
        return "GenomePosition({}, {})".format(*map(str, [self.chrom, self.pos]))

    def __repr__(self):
        return str(self)


class OverlappingWindowSweeper:
    """Sweeping window helper with overlapping windows

    Given a window length ``w`` and an overlap length ``o``, collects records until at least the
    window ``w`` is filled up.  Then, for all but the first and last window on a chromosome,
    records in in the inner ``w - 2*o`` are processed (for the first window, the first ``o`` and
    for the last window, the last ``o`` are processed as well).

    Windows falling out of the window, are yielded from the ``push()`` function.

    As soon as a window is full, each record is processed with the other records in its
    environment.  The call will be ``callback(record, left, right)`` where ``left`` and ``right``
    are iterables with the records that lie left and right of ``record``.

    Note that the callback is only performed once for each record.  The important guarantee here is
    that the callback is called for each record with all records left and right of it in the
    overlap region.

    You have to call ``finish()`` to process the remaining records in the buffer.
    """

    def __init__(self, window_length, overlap, callback, converter=None):
        #: Window length to use, in base pairs
        self.window_length = window_length
        #: Overlap to use, in base pairs
        self.overlap = overlap
        #: Callback to use for processing records
        self.callback = callback
        #: Optional converter to GenomePosition
        self.converter = converter or (lambda x: x)
        #: Current buffer of records
        self.buf = []
        #: Positions for the buffer of records
        self.pos = []
        self._start_chrom(None)

    def _start_chrom(self, name):
        self.chrom = name
        self.begin_outer = 0
        self.begin_inner = 0
        self.end_inner = self.window_length
        self.end_outer = self.window_length + self.overlap

    def push(self, obj, converter=None):
        """Push an object into the window with an optional explicit converter

        Objects falling out of the sweeper will be yielded.
        """
        converter = converter or self.converter
        pos = converter(obj)
        if self.chrom != pos.chrom:
            # Case: first of new chromosome
            yield from self._handle_new_window(pos, obj)
            self._start_chrom(pos.chrom)
            self.pos = [pos]
            self.buf = [obj]
        elif pos >= GenomePosition(self.chrom, self.end_outer):
            # Case: right of outer window
            yield from self._handle_new_window(pos, obj)
            self.pos.append(pos)
            self.buf.append(obj)
        else:
            # Base case: inside of outer window, just append
            self.pos.append(pos)
            self.buf.append(obj)

    def _handle_new_window(self, pos, obj):
        while self.pos and not pos.in_range(self.chrom, self.begin_outer, self.end_outer):
            # Get start position in self.pos and self.obj for processing
            if self.begin_outer == 0:
                # Start at first record
                start = 0
            else:
                # Search for first record in inner window
                query = GenomePosition(self.chrom, self.begin_inner)
                start = bisect.bisect_left(self.pos, query)
            # Get end position in self.pos and self.obj for processing
            if pos.chrom != self.chrom:
                # End at last record
                end = len(self.pos)
            else:
                # Search for end point of inner window
                query = GenomePosition(self.chrom, self.end_inner)
                end = bisect.bisect_right(self.pos, query)
            # Process records inside inner window
            for i in range(start, end):
                self.callback(self.buf[i], self.buf[:i], self.buf[i + 1 :])
            # Advance in buffers
            yield from self.buf[:end]
            self.buf = self.buf[end:]
            self.pos = self.pos[end:]
            # Advance window to next
            self.begin_inner = self.end_inner
            self.begin_outer = max(self.begin_inner - self.overlap, 0)
            self.end_inner += self.window_length
            self.end_outer = self.end_inner + self.overlap

    def finish(self, sentinel="<sentinel>"):
        """Process all remaining records

        ``sentinel`` will be used as the chromosome name of a sentinel position, must not be a
        previously used chromosome identifier.
        """
        yield from self.push(GenomePosition(sentinel, 0), lambda x: x)
