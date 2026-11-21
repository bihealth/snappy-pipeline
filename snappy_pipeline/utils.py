# -*- coding: utf-8 -*-
"""Utility code"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import functools
from typing import Callable, Generator, Iterable, TypeVar

F = TypeVar("F")
T = TypeVar("T")


def listify(gen: Callable[..., Generator[T, None, None]]) -> Callable[..., list[T]]:
    """Decorator that converts a generator into a function which returns a list

    Use it in the case where a generator is easier to write but you want
    to enforce returning a list::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield i
    """

    # functools.wraps keeps the signature visible to inspect, which Snakemake uses to decide
    # which keyword arguments an input or params function accepts.
    @functools.wraps(gen)
    def patched(*args, **kwargs) -> list[T]:
        """Wrapper function"""
        return list(gen(*args, **kwargs))

    return patched


def dictify[**P](gen) -> Callable[P, dict]:
    """Decorator that converts a generator into a function which returns a dict

    Use it in the case where a generator is easier to write but you want
    to enforce returning a dict::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield 'key{}'.format(i), i
    """

    @functools.wraps(gen)
    def patched(*args, **kwargs) -> dict[F, T]:
        """Wrapper function"""
        return dict(gen(*args, **kwargs))

    return patched


def flatten(coll: list[str | list[str] | Iterable[str]]) -> list[str]:
    """Flatten collection of strings or list of strings.

    Source: https://stackoverflow.com/a/17865033
    """
    for i in coll:
        if isinstance(i, Iterable) and not isinstance(i, str):
            for subc in flatten(i):
                yield subc
        else:
            yield i
