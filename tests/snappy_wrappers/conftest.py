from collections import namedtuple
from unittest.mock import MagicMock

from pyfakefs import fake_filesystem
import pytest


@pytest.fixture
def fai_file_content():
    """Returns FAI file content based on hs37d5 (chromosome 1 only)."""
    return "1\t249250621\t52\t60\t61"


@pytest.fixture
def fake_fs():
    """Return ``namedtuple`` with fake file system objects."""
    klass = namedtuple("FakeFsBundle", "fs os open inter_process_lock")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    fake_lock = MagicMock()
    return klass(fs=fake_fs, os=fake_os, open=fake_open, inter_process_lock=fake_lock)


@pytest.fixture
def somatic_variant_fake_fs(fake_fs, fai_file_content):
    """Return fake file system setup with files for the somatic variant calling workflow."""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create static files
    fake_fs.fs.create_file("/path/to/ref.fa", create_missing_dirs=True)
    fake_fs.fs.create_file("/path/to/ref.fa.fai", contents=fai_file_content)
    return fake_fs


def patch_module_fs(module_name, fake_fs, mocker):
    """Helper function to mock out the file-system related things in the module with the given
    name using the given fake_fs and pytest-mock mocker
    """
    mocker.patch("{}.os".format(module_name), fake_fs.os)
    mocker.patch("{}.open".format(module_name), fake_fs.open, create=True)
    mocker.patch("{}.os".format(module_name), fake_fs.os)
