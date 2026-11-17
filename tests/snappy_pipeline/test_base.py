# -*- coding: utf-8 -*-
"""Tests for base module code"""

import pytest

from snappy_pipeline.base import (
    ConfigError,
    UnknownFiltrationSourceException,
    UnsupportedActionException,
    in_task,
)
from snappy_pipeline.workflow_model import TaskModel


def test_config_error_names_the_task_it_is_attached_to():
    """The message states the problem; ``in_task`` attaches the task once, innermost wins."""
    mapping = TaskModel(name="mapping", step="ngs_mapping")
    calling = TaskModel(name="calling", step="variant_calling")
    assert str(ConfigError("no reads")) == "no reads"
    with pytest.raises(ConfigError) as error:
        with in_task(calling):
            with in_task(mapping):
                raise ConfigError("no reads")
    assert error.value.task is mapping
    assert str(error.value) == "Task 'mapping' (step ngs_mapping): no reads"
    assert error.value.message == "no reads"


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
