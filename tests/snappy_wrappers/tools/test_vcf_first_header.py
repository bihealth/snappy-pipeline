# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.vcf_first_header``"""

from io import StringIO
import textwrap
from unittest.mock import patch

from snappy_wrappers.tools import vcf_first_header


def test_vcf_first_header(capsys):
    fake_stdin = StringIO(
        textwrap.dedent(
            r"""
    #CHROM
    line1
    #CHROM
    line2
    """
        ).lstrip()
    )

    with patch("snappy_wrappers.tools.vcf_first_header.sys.stdin", fake_stdin):
        vcf_first_header.main()

    out, err = capsys.readouterr()
    assert (
        out
        == textwrap.dedent(
            r"""
    #CHROM
    line1
    line2
    """
        ).lstrip()
    )
    assert err == ""
