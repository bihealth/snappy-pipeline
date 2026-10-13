# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.vcf_first_header``"""

import io
import textwrap

import pytest

from snappy_wrappers.tools import genome_windows


@pytest.fixture
def ref_fasta_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake BED files
    fcontents = textwrap.dedent(
        """
    1\t249250621\t52\t60\t61
    2\t243199373\t253404903\t60\t61
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/ref.fasta.fai", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_genome_windows(capsys, ref_fasta_fake_fs, mocker):
    mocker.patch("argparse.open", ref_fasta_fake_fs.open, create=True)
    mocker.patch("argparse._os", ref_fasta_fake_fs.os)

    genome_windows.main(
        ["--fai-file", "/work/ref.fasta.fai", "--window-size", str(100 * 1000 * 1000)]
    )

    out, err = capsys.readouterr()
    assert (
        out
        == textwrap.dedent(
            r"""
    1:1-100,000,000
    1:100,000,001-200,000,000
    1:200,000,001-249,250,621
    2:1-100,000,000
    2:100,000,001-200,000,000
    2:200,000,001-243,199,373
    """
        ).lstrip()
    )
    assert err == ""


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
