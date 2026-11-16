# -*- coding: utf-8 -*-
"""Tests for ``yaml_utils``."""

from snappy_pipeline.apps.impl.yaml_utils import remove_yaml_comment_lines


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
