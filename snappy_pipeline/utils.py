# -*- coding: utf-8 -*-
"""Utility code"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

from typing import Any, Callable, Generator, Iterable, TypeVar

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


def try_or_none(func, exceptions):
    """Helper that tries to execute the function

    If one of the ``exceptions`` is raised then return None
    """
    try:
        return func()
    except exceptions:
        return None


def is_none(value: Any) -> bool:
    """Helper function returning whether ``value is None``"""
    return value is None


def is_not_none(value: Any) -> bool:
    """Helper function returning whether ``value is not None``"""
    return value is not None


class DictQuery(dict):
    """Helper class for comfortable access to nested dicts with ``str`` keys.

    Source:

    - https://www.haykranen.nl/2016/02/13/handling-complex-nested-dicts-in-python/
    """

    def get(self, path: str, default=None):
        keys = path.split("/")
        val = None

        for key in keys:
            if val:
                if isinstance(val, list):
                    val = [v.get(key, default) if v else None for v in val]
                else:
                    val = val.get(key, default)
            else:
                val = dict.get(self, key, default)

            if not val:
                break

        return val
