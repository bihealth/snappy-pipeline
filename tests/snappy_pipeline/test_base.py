# -*- coding: utf-8 -*-
"""Tests for base module code"""

import pytest

from snappy_pipeline.base import (
    InvalidConfiguration,
    MissingConfiguration,
    UnknownFiltrationSourceException,
    UnsupportedActionException,
)


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
