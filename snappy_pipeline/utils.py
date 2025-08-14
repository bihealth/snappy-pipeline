# -*- coding: utf-8 -*-
"""Utility code"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import typing


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
        """Wrapper function"""
        return list(gen(*args, **kwargs))

    return patched


def dictify[**P](gen) -> typing.Callable[P, dict]:
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
        """Wrapper function"""
        return dict(gen(*args, **kwargs))

    return patched


def flatten(coll: typing.List[typing.Union[str, typing.List[str]]]) -> typing.List[str]:
    """Flatten collection of strings or list of strings.

    Source: https://stackoverflow.com/a/17865033
    """
    for i in coll:
        if isinstance(i, typing.Iterable) and not isinstance(i, str):
            for subc in flatten(i):
                yield subc
        else:
            yield i


def try_or_none(func, exceptions):
    """Helper that tries to execute the function

    If one of the ``exceptions`` is raised then return None
    """
    try:
        return func()
    except exceptions:
        return None


def is_none(value):
    """Helper function returning whether ``value is None``"""
    return value is None


def is_not_none(value):
    """Helper function returning whether ``value is not None``"""
    return value is not None

