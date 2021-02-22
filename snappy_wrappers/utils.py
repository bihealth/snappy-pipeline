# -*- coding: utf-8 -*-
"""Utility code for snappy_wrappers"""

import os

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def basename(fname, strip_suffix=""):
    """Unix `basename`-like function"""
    fname = os.path.basename(fname)
    if strip_suffix and fname.endswith(strip_suffix):
        return fname[: -len(strip_suffix)]
    else:
        return fname


def listify(gen):
    """Decorator that converts a generator into a function which returns a list

    Use it in the case where a generator is easier to write but you want
    to enforce returning a list::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield i
    """

    def patched(*args, **kwargs):
        """Wrapper function """
        return list(gen(*args, **kwargs))

    return patched


def dictify(gen):
    """Decorator that converts a generator into a function which returns a dict

    Use it in the case where a generator is easier to write but you want
    to enforce returning a dict::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield 'key{}'.format(i), i
    """

    def patched(*args, **kwargs):
        """Wrapper function """
        return dict(gen(*args, **kwargs))

    return patched
