# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.bed_filter_jaccard``"""

import pytest

from snappy_wrappers.tools import bed_filter_jaccard


@pytest.fixture
def input_txt_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = "1\t1\t100\t1\t10\t110\t90\n2\t1\t200\t.\t-1\t-1\t0\n"

    fake_fs.fs.create_file("/work/input.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_bed_filter_jaccard_subtract(capsys, input_txt_fake_fs, mocker):
    mocker.patch("argparse.open", input_txt_fake_fs.open, create=True)
    mocker.patch("argparse._os", input_txt_fake_fs.os)

    bed_filter_jaccard.main(["--input-file", "/work/input.txt", "--operation", "subtract"])

    expected_out = "2\t1\t200\t.\t-1\t-1\t0\n"

    out, err = capsys.readouterr()
    assert out == expected_out
    assert err == ""


def test_bed_filter_jaccard_intersect(capsys, input_txt_fake_fs, mocker):
    mocker.patch("argparse.open", input_txt_fake_fs.open, create=True)
    mocker.patch("argparse._os", input_txt_fake_fs.os)

    bed_filter_jaccard.main(["--input-file", "/work/input.txt", "--operation", "intersect"])

    expected_out = "1\t1\t100\t1\t10\t110\t90\n"

    out, err = capsys.readouterr()
    assert out == expected_out
    assert err == ""
