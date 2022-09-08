# -*- coding: utf-8 -*-
"""Tests for ``yaml_utils``."""
import ruamel.yaml as ruamel_yaml

import pytest

from snappy_pipeline.apps.impl.yaml_utils import remove_non_required, remove_yaml_comment_lines


@pytest.fixture
def yaml_obj():
    return ruamel_yaml.YAML()


@pytest.fixture
def step_config(yaml_obj):
    """Returns genetic step config"""
    return yaml_obj.load(
        remove_yaml_comment_lines(
            r"""
            step_config:
              step_name:
                tool:
                  argument1: auto
                  argument2: 16
                  argument3: 8
    """
        )
    )


@pytest.fixture
def step_config_required_tag(yaml_obj):
    """Returns genetic step config"""
    return yaml_obj.load(
        remove_yaml_comment_lines(
            r"""
            step_config:
              step_name:
                tool:
                  argument1: REQUIRED # Shouldn't work if not in the comments
                  argument2: null
                  argument3: auto
                  argument4: 16       # Required in the comments
                  argument5: 8
    """
        )
    )


@pytest.fixture
def step_config_optional_tag(yaml_obj):
    """Returns genetic step config"""
    return yaml_obj.load(
        remove_yaml_comment_lines(
            r"""
            step_config:
              step_name:
                tool:
                  argument1: null
                  argument2: auto  # Optional in the comments
                  argument3: 16
                  argument4: 8
    """
        )
    )


@pytest.fixture
def step_config_both_tags(yaml_obj):
    """Returns genetic step config"""
    return yaml_obj.load(
        remove_yaml_comment_lines(
            r"""
            step_config:
              step_name:
                tool:
                  argument1: null
                  argument2: null
                  argument3: REQUIRED  # Shouldn't work if not in the comments
                  argument4: 16        # REQUIRED in the comments
                  argument5: 8         # OPTIONAL in the comments
    """
        )
    )


def test_remove_non_required_no_tags(step_config):
    """Tests remove_non_required() - config has no tags"""
    actual = remove_non_required(yaml_obj=step_config)
    assert len(actual) == 0


def test_remove_non_required_required_tag(step_config_required_tag):
    """Tests remove_non_required() - config has REQUIRED tags"""
    expected = {"argument4": 16}
    result = remove_non_required(yaml_obj=step_config_required_tag)
    observed = result["step_config"]["step_name"]["tool"]
    assert observed == expected


def test_remove_non_required_optional_tag(step_config_optional_tag):
    """Tests remove_non_required() - config has REQUIRED tags"""
    expected = {"argument2": "auto"}
    result = remove_non_required(yaml_obj=step_config_optional_tag)
    observed = result["step_config"]["step_name"]["tool"]
    assert observed == expected


def test_remove_non_required_both_tags(step_config_both_tags):
    """Tests remove_non_required() - config has both REQUIRED and OPTIONAL tags"""
    expected = {"argument4": 16, "argument5": 8}
    result = remove_non_required(yaml_obj=step_config_both_tags)
    observed = result["step_config"]["step_name"]["tool"]
    assert observed == expected


def test_remove_yaml_comment_lines():
    config_str = r"""
            step_config:
              # Commented line
              step_name:
                tool:
                  # Commented line
                  argument1: auto
                  argument2: 16
                  #
                  argument3: 8
    """
    actual = remove_yaml_comment_lines(yaml_str=config_str)
    assert "#" not in actual
