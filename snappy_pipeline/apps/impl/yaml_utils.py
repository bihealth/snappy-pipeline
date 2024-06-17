# -*- coding: utf-8 -*-
"""Utilities for processing YAML configuration"""

from collections.abc import MutableMapping, MutableSequence
import re

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


def remove_non_required(yaml_obj):
    """Remove non required arguments

    Remove items that are not marked as 'required' nor 'optional' with a comment (case-insensitive).

    :param yaml_obj: YAML configuration.
    :type yaml_obj: ruamel.yaml.comments.CommentedMap

    :return: Returns ordered dictionary-like with arguments to be included in the snappy
    directory config. Output is either ``CommentedMap`` or ``CommentedSeq``.
    """
    if isinstance(yaml_obj, (dict, MutableMapping)):
        result = CommentedMap()
        for key, value in yaml_obj.items():
            required = key in yaml_obj.ca.items and (
                "required" in yaml_obj.ca.items[key][2].value.lower()
                or "optional" in yaml_obj.ca.items[key][2].value.lower()
            )

            if isinstance(value, (dict, MutableMapping, list, MutableSequence)):
                value = remove_non_required(value)
                required = required or value
            if required:
                result[key] = value
                if key in yaml_obj.ca.items:
                    result.yaml_add_eol_comment(yaml_obj.ca.items[key][2].value, key)
        return result
    elif isinstance(yaml_obj, (list, MutableSequence)):
        result = CommentedSeq()
        for key, value in enumerate(yaml_obj):
            comment = ""
            if yaml_obj.ca.items and key in yaml_obj.ca.items and yaml_obj.ca.items[key][2]:
                comment = yaml_obj.ca.items[key][2].value
            required = key in yaml_obj.ca.items and (
                "required" in comment.lower() or "optional" in comment.lower()
            )
            if isinstance(value, (dict, MutableMapping, list, MutableSequence)):
                value = remove_non_required(value)
                required = required or value
            if required:
                result.append(value)
                if key in yaml_obj.ca.items:
                    result.yaml_add_eol_comment(yaml_obj.ca.items[key][2].value)
        return result
    else:
        assert False, "Input must be either dict or list."
