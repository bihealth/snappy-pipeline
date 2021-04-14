# -*- coding: utf-8 -*-
"""Tests for base module code"""

from snappy_pipeline.base import merge_kwargs


def test_merge_kwargs():
    """Tests dictionary merger for shortcut sheet keyword arguments."""
    # Initialise variables
    global_kwargs = {1: "one", 2: "two"}
    sheet_kwargs = {3: "three"}
    merged_kwargs = {1: "one", 2: "two", 3: "three"}

    # None if both none
    expected = None
    actual = merge_kwargs(first_kwargs=None, second_kwargs=None)
    assert actual == expected

    # Only sheet if global is none
    actual = merge_kwargs(first_kwargs=None, second_kwargs=sheet_kwargs)
    assert actual == sheet_kwargs

    # Only global if sheet is none
    actual = merge_kwargs(first_kwargs=global_kwargs, second_kwargs=None)
    assert actual == global_kwargs

    # Merged if both presents
    actual = merge_kwargs(first_kwargs=global_kwargs, second_kwargs=sheet_kwargs)
    assert actual == merged_kwargs
