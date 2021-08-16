# -*- coding: utf-8 -*-
"""Basic utility code for snappy_pipeline
"""

from collections import OrderedDict
from collections.abc import MutableMapping
from copy import deepcopy
import os
import sys
import warnings

import ruamel.yaml as yaml

# TODO: This has to go away once biomedsheets is a proper, halfway-stable module
try:
    from biomedsheets.ref_resolver import RefResolver
except ImportError:
    warnings.warn("module biomedsheets not found", UserWarning)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class InvalidConfiguration(Exception):
    """Raised on invalid configuration"""


class MissingConfiguration(InvalidConfiguration):
    """Raised on missing configuration"""


class UnsupportedActionException(Exception):
    """Raised when user try to call action that isn't supported."""


class UnknownFiltrationSourceException(Exception):
    """Raised when user try to request an unknown filtration source."""


def expand_ref(config_path, dict_data, lookup_paths=None, dict_class=OrderedDict):
    """Expand "$ref" in JSON-like data ``dict_data``

    Returns triple:

    - path to resolved file
    - paths containing included config files
    - config files included
    """
    lookup_paths = lookup_paths or [os.getcwd()]
    resolver = RefResolver(lookup_paths=lookup_paths, dict_class=dict_class)
    # Perform resolution
    resolved = resolver.resolve("file://" + config_path, dict_data)
    # Collect paths of all included configuration files, important for
    # data set importing later on
    lookup_paths = list(lookup_paths)  # copy!
    config_files = []  # config files (not URLs) read
    for url in resolver.cache:
        if url.startswith("file://"):
            config_files.append(os.path.abspath(url[len("file://") :]))
            dirname = os.path.dirname(url[len("file://") :])
            if not dirname:
                dirname = "."
            if dirname not in lookup_paths:
                lookup_paths.append(dirname)
    return resolved, tuple(lookup_paths), tuple(config_files)


def _ordered_dump(data, stream=None, dumper_klass=yaml.Dumper, **kwargs):
    """Helper function for ordered dumping of structs as YAML"""

    class OrderedDumper(dumper_klass):
        """Helper class"""

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items()
        )

    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwargs)


def print_config(config, file=sys.stderr):
    """Print human-readable version of configuration to ``file``"""
    print("\nConfiguration", file=file)
    print("-------------\n", file=file)
    _ordered_dump(config, stream=file, default_flow_style=False)


def print_sample_sheets(step, file=sys.stderr):
    """Print loaded sample sheets from ``BaseStep`` in human-readable format"""
    for info in step.data_set_infos:
        print("\nSample Sheet {}".format(info.sheet_path), file=file)
        print("-------------" + "-" * len(info.sheet_path) + "\n", file=file)
        _ordered_dump(info.sheet.json_data, stream=file, default_flow_style=False)


def merge_kwargs(first_kwargs, second_kwargs):
    """Merge two keyword arguments.

    :param first_kwargs: First keyword arguments dictionary.
    :type first_kwargs: dict

    :param second_kwargs: Second keyword arguments dictionary.
    :type second_kwargs: dict

    :return: Returns merged dictionary with inputted keyword arguments.
    """
    # Global if no individual dict
    if first_kwargs and (not second_kwargs):
        return first_kwargs
    # Individual if no global dict
    elif (not first_kwargs) and second_kwargs:
        return second_kwargs
    # Merge dicts if both defined
    elif first_kwargs and second_kwargs:
        global_copy_kwargs = deepcopy(first_kwargs)
        global_copy_kwargs.update(second_kwargs)
        return global_copy_kwargs
    # None if both None
    else:
        return None


def merge_dicts(dict1, dict2, dict_class=OrderedDict):
    """Merge dictionary ``dict2`` into ``dict1``"""

    def _merge_inner(dict1, dict2):
        for k in set(dict1.keys()).union(dict2.keys()):
            if k in dict1 and k in dict2:
                if isinstance(dict1[k], (dict, MutableMapping)) and isinstance(
                    dict2[k], (dict, MutableMapping)
                ):
                    yield k, dict_class(_merge_inner(dict1[k], dict2[k]))
                else:
                    # If one of the values is not a dict, you can't continue
                    # merging it.  Value from second dict overrides one in
                    # first and we move on.
                    yield k, dict2[k]
            elif k in dict1:
                yield k, dict1[k]
            else:
                yield k, dict2[k]

    return dict_class(_merge_inner(dict1, dict2))


def snakefile_path(step_name):
    """Return absolute path to Snakefile for the given step name"""
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "workflows", step_name, "Snakefile")
    )
