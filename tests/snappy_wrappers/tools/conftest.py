# -*- coding: utf-8 -*-
"""Shared fixtures for the unit tests"""

from collections import namedtuple

from pyfakefs import fake_filesystem
import pytest


@pytest.fixture
def fake_fs():
    """Return ``namedtuple`` with fake file system objects"""
    klass = namedtuple("FakeFsBundle", "fs os open")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    return klass(fs=fake_fs, os=fake_os, open=fake_open)
