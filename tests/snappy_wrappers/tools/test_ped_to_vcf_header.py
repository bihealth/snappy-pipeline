# -*- coding: utf-8 -*-

import os

from snappy_wrappers.tools import ped_to_vcf_header

import pytest
import textwrap
import tempfile
from collections import namedtuple


Donor = namedtuple("Donor", ["family", "id", "father", "mother", "sex", "disease"])

DONORS_COMPLETE = [
    Donor(*["FAM", "II-1", "I-1", "I-2", "1", "2"]),
    Donor(*["FAM", "I-1", "0", "0", "1", "1"]),
    Donor(*["FAM", "I-2", "0", "0", "2", "1"]),
]

EXPECTED_COMPLETE = textwrap.dedent(
    r"""
##META=<ID=Sex,Type=String,Number=1,Values=[Unknown, Male, Female]>
##META=<ID=Disease,Type=String,Number=1,Values=[Unknown, Unaffected, Affected]>
##SAMPLE=<ID=II-1,Sex=Male,Disease=Affected>
##SAMPLE=<ID=I-1,Sex=Male,Disease=Unaffected>
##SAMPLE=<ID=I-2,Sex=Female,Disease=Unaffected>
##PEDIGREE=<ID=II-1,Family=FAM,Father=I-1,Mother=I-2>
"""
).strip()


@pytest.fixture
def ped_file_complete(fake_fs):
    """Return fake file system with PED file"""

    content = textwrap.dedent(
        """
    # comment
    FAM	II-1\tI-1\tI-2\t1\t2
    FAM I-1     0\t0\t1\t1
    FAM I-2     0\t0\t2\t1
    """
    ).strip()

    fake_fs.fs.create_file("/test.ped", create_missing_dirs=True, contents=content)

    return fake_fs


def test_parse_ped_complete(ped_file_complete, mocker):
    mocker.patch("%s.open" % __name__, ped_file_complete.open, create=True)

    with open("/test.ped", "rt") as ped:
        assert DONORS_COMPLETE == list(ped_to_vcf_header.parse_ped(ped))


def test_ped_vcf_header_complete():
    output = ped_to_vcf_header.ped_vcf_header(DONORS_COMPLETE)

    assert EXPECTED_COMPLETE == output


def test_write_header_snippet_complete():
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")
        ped_to_vcf_header.write_header_snippet(EXPECTED_COMPLETE, output_path)

        with open(output_path, "r") as fh:
            assert EXPECTED_COMPLETE == fh.read().rstrip()


def test_main_complete(ped_file_complete, mocker):
    mocker.patch("argparse.open", ped_file_complete.open, create=True)
    mocker.patch("argparse._os", ped_file_complete.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")

        ped_to_vcf_header.main(["--ped-file", "/test.ped", "--output", output_path])

        assert os.path.exists(output_path)


DONORS_PARENT_MISSING = [
    Donor(*["FAM", "II-1", "I-1", "0", "1", "2"]),
    Donor(*["FAM", "I-1", "0", "0", "1", "1"]),
]

EXPECTED_PARENT_MISSING = textwrap.dedent(
    r"""
##META=<ID=Sex,Type=String,Number=1,Values=[Unknown, Male, Female]>
##META=<ID=Disease,Type=String,Number=1,Values=[Unknown, Unaffected, Affected]>
##SAMPLE=<ID=II-1,Sex=Male,Disease=Affected>
##SAMPLE=<ID=I-1,Sex=Male,Disease=Unaffected>
##PEDIGREE=<ID=II-1,Family=FAM,Father=I-1,Mother=0>
"""
).strip()


@pytest.fixture
def ped_file_parent_missing(fake_fs):
    """Return fake file system with PED file"""

    content = textwrap.dedent(
        """
    # comment
    FAM	II-1\tI-1\t0\t1\t2
    FAM I-1     0\t0\t1\t1
    """
    ).strip()

    fake_fs.fs.create_file("/test.ped", create_missing_dirs=True, contents=content)

    return fake_fs


def test_parse_ped_parent_missing(ped_file_parent_missing, mocker):
    mocker.patch("%s.open" % __name__, ped_file_parent_missing.open, create=True)

    with open("/test.ped", "rt") as ped:
        assert DONORS_PARENT_MISSING == list(ped_to_vcf_header.parse_ped(ped))


def test_ped_vcf_header_parent_missing():
    output = ped_to_vcf_header.ped_vcf_header(DONORS_PARENT_MISSING)

    assert EXPECTED_PARENT_MISSING == output


def test_write_header_snippet_parent_missing():
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")
        ped_to_vcf_header.write_header_snippet(EXPECTED_PARENT_MISSING, output_path)

        with open(output_path, "r") as fh:
            assert EXPECTED_PARENT_MISSING == fh.read().rstrip()


def test_main_parent_missing(ped_file_parent_missing, mocker):
    mocker.patch("argparse.open", ped_file_parent_missing.open, create=True)
    mocker.patch("argparse._os", ped_file_parent_missing.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")

        ped_to_vcf_header.main(["--ped-file", "/test.ped", "--output", output_path])

        assert os.path.exists(output_path)


DONORS_SINGLETON = [Donor(*["FAM", "II-1", "0", "0", "1", "2"])]

EXPECTED_SINGLETON = textwrap.dedent(
    r"""
##META=<ID=Sex,Type=String,Number=1,Values=[Unknown, Male, Female]>
##META=<ID=Disease,Type=String,Number=1,Values=[Unknown, Unaffected, Affected]>
##SAMPLE=<ID=II-1,Sex=Male,Disease=Affected>
##PEDIGREE=<ID=II-1,Family=FAM,Father=0,Mother=0>
"""
).strip()


@pytest.fixture
def ped_file_singleton(fake_fs):
    """Return fake file system with PED file"""

    content = textwrap.dedent(
        """
    # comment
    FAM	II-1\t0\t0\t1\t2
    """
    ).strip()

    fake_fs.fs.create_file("/test.ped", create_missing_dirs=True, contents=content)

    return fake_fs


def test_parse_ped_singleton(ped_file_singleton, mocker):
    mocker.patch("%s.open" % __name__, ped_file_singleton.open, create=True)

    with open("/test.ped", "rt") as ped:
        assert DONORS_SINGLETON == list(ped_to_vcf_header.parse_ped(ped))


def test_ped_vcf_header_singleton():
    output = ped_to_vcf_header.ped_vcf_header(DONORS_SINGLETON)

    assert EXPECTED_SINGLETON == output


def test_write_header_snippet_singleton():
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")
        ped_to_vcf_header.write_header_snippet(EXPECTED_SINGLETON, output_path)

        with open(output_path, "r") as fh:
            assert EXPECTED_SINGLETON == fh.read().rstrip()


def test_main_singleton(ped_file_singleton, mocker):
    mocker.patch("argparse.open", ped_file_singleton.open, create=True)
    mocker.patch("argparse._os", ped_file_singleton.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "header")

        ped_to_vcf_header.main(["--ped-file", "/test.ped", "--output", output_path])

        assert os.path.exists(output_path)
