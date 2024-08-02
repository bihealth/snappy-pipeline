# -*- coding: utf-8 -*-
"""Tests for ``snappy_wrappers.tools.quickvenn``"""

import os
import tempfile
import textwrap

import pytest

from snappy_wrappers.tools import quickvenn


@pytest.fixture
def one_set_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = textwrap.dedent(
        """
    A\tlabel
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/set.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_quickvenn_one_set(one_set_fake_fs, mocker):
    mocker.patch("argparse.open", one_set_fake_fs.open, create=True)
    mocker.patch("argparse._os", one_set_fake_fs.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "output.png")
        quickvenn.main(["--input-shared-counts", "/work/set.txt", "--output-image", output_path])

        assert os.path.exists(output_path)


@pytest.fixture
def two_sets_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = textwrap.dedent(
        """
    A\tlabel
    B\tlabel
    A\tB\tlabel
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/set.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_quickvenn_two_sets(two_sets_fake_fs, mocker):
    mocker.patch("argparse.open", two_sets_fake_fs.open, create=True)
    mocker.patch("argparse._os", two_sets_fake_fs.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "output.png")
        quickvenn.main(["--input-shared-counts", "/work/set.txt", "--output-image", output_path])

        assert os.path.exists(output_path)


@pytest.fixture
def three_sets_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = textwrap.dedent(
        """
    A\tlabel
    B\tlabel
    A\tB\tlabel
    C\tlabel
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/set.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_quickvenn_three_sets(three_sets_fake_fs, mocker):
    mocker.patch("argparse.open", three_sets_fake_fs.open, create=True)
    mocker.patch("argparse._os", three_sets_fake_fs.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "output.png")
        quickvenn.main(["--input-shared-counts", "/work/set.txt", "--output-image", output_path])

        assert os.path.exists(output_path)


@pytest.fixture
def four_sets_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = textwrap.dedent(
        """
    A\tlabel
    B\tlabel
    A\tB\tlabel
    C\tlabel
    D\tlabel
    C\tD\tlabel
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/set.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_quickvenn_four_sets(four_sets_fake_fs, mocker):
    mocker.patch("argparse.open", four_sets_fake_fs.open, create=True)
    mocker.patch("argparse._os", four_sets_fake_fs.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "output.png")
        quickvenn.main(["--input-shared-counts", "/work/set.txt", "--output-image", output_path])

        assert os.path.exists(output_path)


@pytest.fixture
def five_sets_fake_fs(fake_fs):
    """Return fake file system setup with files with an FAI file"""
    # Create fake bedtools intersect -wao output
    fcontents = textwrap.dedent(
        """
    A\tlabel
    B\tlabel
    A\tB\tlabel
    C\tlabel
    D\tlabel
    C\tD\tlabel
    E\tlabel
    """
    ).lstrip()

    fake_fs.fs.create_file("/work/set.txt", create_missing_dirs=True, contents=fcontents)
    return fake_fs


def test_quickvenn_five_sets(five_sets_fake_fs, mocker):
    mocker.patch("argparse.open", five_sets_fake_fs.open, create=True)
    mocker.patch("argparse._os", five_sets_fake_fs.os)

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_path = os.path.join(tmpdirname, "output.png")
        quickvenn.main(["--input-shared-counts", "/work/set.txt", "--output-image", output_path])

        assert os.path.exists(output_path)
