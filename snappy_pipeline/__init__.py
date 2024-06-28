# -*- coding: utf-8 -*-

from .base import expand_ref, merge_dictlikes, print_config, print_sample_sheets

__author__ = """Manuel Holtgrewe"""
__email__ = "manuel.holtgrewe@bih-charite.de"

from snappy_pipeline._version import __version__

__all__ = ["__version__", "expand_ref", "merge_dicts", "print_config", "print_sample_sheets"]
