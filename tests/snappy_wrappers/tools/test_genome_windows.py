# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.genome_windows``"""

import io

from snappy_wrappers.tools import genome_windows


def test_matches_any_uses_regular_expressions_with_search_semantics():
    patterns = ["^hs37d5$", "_decoy$", "^HLA-"]
    assert genome_windows.matches_any("hs37d5", patterns)
    assert genome_windows.matches_any("chr1_KI270_decoy", patterns)
    assert genome_windows.matches_any("HLA-A*01:01", patterns)
    assert not genome_windows.matches_any("hs37d5x", patterns)
    assert not genome_windows.matches_any("chr1", patterns)
    assert not genome_windows.matches_any("my_HLA-A", patterns)
    assert not genome_windows.matches_any("chr1", [])


def test_yield_contigs_skips_ignored_contigs():
    fai = io.StringIO("1\t10\t0\t60\t61\nhs37d5\t10\t0\t60\t61\nGL000220.1\t10\t0\t60\t61\n")
    assert list(genome_windows.yield_contigs(fai, ["^hs37d5$", r"^GL\d+"])) == ["1"]
