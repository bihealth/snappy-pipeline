# -*- coding: utf-8 -*-
"""Utilities for processing YAML configuration"""

import re


def remove_yaml_comment_lines(yaml_str):
    """
    :param yaml_str: YAML configuration as string. Expected use case: clean DEFAULT_CONFIG from
    workflows.
    :type yaml_str: str

    :return: Returns default configuration YAML string without commented nor empty lines.
    """
    result = []
    for line in yaml_str.splitlines(True):
        if not (re.match(r"^\s*#", line) or line == "\n"):
            result.append(line)
    return "".join(result)
