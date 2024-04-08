# -*- coding: utf-8 -*-

from .base import expand_ref, merge_dicts, print_config, print_sample_sheets

# from ._version import get_versions

__author__ = """Manuel Holtgrewe"""
__email__ = "manuel.holtgrewe@bih-charite.de"
__version__ = "master"

# del get_versions

from . import _version

__version__ = _version.get_versions()["version"]
__all__ = ["expand_ref", "merge_dicts", "print_config", "print_sample_sheets"]
