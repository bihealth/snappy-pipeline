# -*- coding: utf-8 -*-
"""Utilities for processing YAML configuration
"""

from collections.abc import MutableSequence, MutableMapping
import re

from ruamel.yaml.comments import CommentedMap, CommentedSeq


def remove_yaml_comment_lines(yaml_str):
    """Get default configuration YAML string without comment-online lines"""
    result = []
    for line in yaml_str.splitlines(True):
        if not re.match(r"^\s*#", line):
            result.append(line)
    return "".join(result)


def remove_non_required(yaml_obj):
    """Remove items that are not marked as "required" with a comment (case-insensitive)"""
    if isinstance(yaml_obj, (dict, MutableMapping)):
        result = CommentedMap()
        for key, value in yaml_obj.items():
            required = (
                key in yaml_obj.ca.items and "required" in yaml_obj.ca.items[key][2].value.lower()
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
            required = key in yaml_obj.ca.items and "required" in comment.lower()
            if isinstance(value, (dict, MutableMapping, list, MutableSequence)):
                value = remove_non_required(value)
                required = required or value
            if required:
                result.append(value)
                if key in yaml_obj.ca.items:
                    result.yaml_add_eol_comment(yaml_obj.ca.items[key][2].value)
        return result
    else:
        assert False, "Must be dict/list"
