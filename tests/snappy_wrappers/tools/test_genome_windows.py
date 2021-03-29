# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.vcf_first_header``"""

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
