# -*- coding: utf-8 -*-
"""Basic support code for the snappy_wrappers"""

import argparse
import inspect
import os

from .utils import dictify

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class WrapperInfo:
    """Base class for generating relevant entries that are passed to the sequential wrapper"""

    def __init__(self, args):
        #: parsed command line arguments
        self.args = args

    @dictify
    def get_value_prefixes(self):
        """Return dict with directive => value prefix

        This can be used for prefixing values with '**' or 'unpack'.  You don't need to override
        this method.
        """
        for key in "input", "output", "params", "resources":
            if isinstance(getattr(self, "get_" + key)(), dict):
                yield key, "**"
            else:
                yield key, ""

    def get_input(self):  # NOSONAR
        """Override to determine what to use for input files"""
        return []

    def get_output(self):  # NOSONAR
        """Override to determine what to use for input files"""
        return []

    def get_params(self):
        """Override to determine what to use for input files

        Make sure to call super class' ``get_params()`` in sub classes and include the keys from
        this class.
        """
        return {"__wrapper": self._get_wrapper()}

    def get_threads(self):  # NOSONAR
        """Return number of threads to use"""
        return 1

    def get_resources(self):  # NOSONAR
        """Return number of threads to use"""
        return {}

    def get_log(self):
        """Return path to log file"""
        if self.args:
            return self.args.log_file or []
        else:
            return []

    def get_config(self):  # NOSONAR
        """Return configuration to use"""
        return {}

    def _get_wrapper(self):
        """Return URL to the wrapper"""
        return "file://" + os.path.dirname(os.path.abspath(inspect.getfile(type(self))))


def build_parser(*args, **kwargs):
    """Build and return basic ArgumentParser"""
    parser = argparse.ArgumentParser(*args, **kwargs)

    parser.add_argument("--verbose", action="store_true", default=False)
    parser.add_argument("--log-file", help="Path to log file")

    return parser
