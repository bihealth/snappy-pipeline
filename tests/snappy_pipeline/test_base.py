# -*- coding: utf-8 -*-
"""Tests for base module code"""

import pytest

from snappy_pipeline.base import (
    InvalidConfiguration,
    MissingConfiguration,
    UnknownFiltrationSourceException,
    UnsupportedActionException,
    merge_kwargs,
)


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


def test_invalid_configuration_exception():
    """Tests InvalidConfiguration raise."""
    error_msg = "Raised InvalidConfiguration"
    with pytest.raises(Exception) as exec_info:
        raise InvalidConfiguration(error_msg)
    assert exec_info.value.args[0] == error_msg


def test_missing_configuration_exception():
    """Tests MissingConfiguration raise."""
    error_msg = "Raised MissingConfiguration"
    with pytest.raises(Exception) as exec_info:
        raise MissingConfiguration(error_msg)
    assert exec_info.value.args[0] == error_msg


def test_unsupported_action_exception():
    """Tests UnsupportedActionException raise."""
    error_msg = "Raised UnsupportedActionException"
    with pytest.raises(Exception) as exec_info:
        raise UnsupportedActionException(error_msg)
    assert exec_info.value.args[0] == error_msg


def test_unknown_filtration_source_exception():
    """Tests UnknownFiltrationSourceException raise."""
    error_msg = "Raised UnknownFiltrationSourceException"
    with pytest.raises(Exception) as exec_info:
        raise UnknownFiltrationSourceException(error_msg)
    assert exec_info.value.args[0] == error_msg
