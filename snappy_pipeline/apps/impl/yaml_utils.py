# -*- coding: utf-8 -*-
"""Utilities for processing YAML configuration"""

import re
from collections.abc import MutableMapping, MutableSequence

from ruamel.yaml.comments import CommentedMap, CommentedSeq


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


def _get_eol_comment_token(yaml_obj, key):
    if hasattr(yaml_obj, "ca") and isinstance(yaml_obj.ca.items, dict):
        item = yaml_obj.ca.items.get(key)
        if item and len(item) > 2:
            return item[2]
    return None


def remove_non_required(yaml_obj):
    """Remove non-required arguments.

    Remove items that are not marked as 'required' nor 'optional' with a comment (case-insensitive).
    """
    if isinstance(yaml_obj, (dict, MutableMapping)):
        result = CommentedMap()
        for key, value in yaml_obj.items():
            comment_token = _get_eol_comment_token(yaml_obj, key)
            comment_str = (
                comment_token.value.lower()
                if (comment_token and hasattr(comment_token, "value") and comment_token.value)
                else ""
            )
            required = "required" in comment_str or "optional" in comment_str

            if isinstance(value, (dict, MutableMapping, list, MutableSequence)):
                value = remove_non_required(value)
                required = required or bool(value)
            if required:
                result[key] = value
                if comment_token and hasattr(comment_token, "value") and comment_token.value:
                    lines = comment_token.value.splitlines()
                    if lines:
                        eol_text = lines[0].strip()
                        if eol_text.startswith("#"):
                            eol_text = eol_text[1:].strip()
                        if eol_text:
                            result.yaml_add_eol_comment(eol_text, key)
        return result
    elif isinstance(yaml_obj, (list, MutableSequence)):
        result = CommentedSeq()
        for key, value in enumerate(yaml_obj):
            comment_token = _get_eol_comment_token(yaml_obj, key)
            comment_str = (
                comment_token.value.lower()
                if (comment_token and hasattr(comment_token, "value") and comment_token.value)
                else ""
            )
            required = "required" in comment_str or "optional" in comment_str
            if isinstance(value, (dict, MutableMapping, list, MutableSequence)):
                value = remove_non_required(value)
                required = required or bool(value)
            if required:
                result.append(value)
                if comment_token and hasattr(comment_token, "value") and comment_token.value:
                    lines = comment_token.value.splitlines()
                    if lines:
                        eol_text = lines[0].strip()
                        if eol_text.startswith("#"):
                            eol_text = eol_text[1:].strip()
                        if eol_text:
                            result.yaml_add_eol_comment(eol_text, key=len(result) - 1)
        return result
    else:
        assert False, "Input must be either dict or list."
