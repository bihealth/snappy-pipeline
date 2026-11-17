# -*- coding: utf-8 -*-
"""Basic utility code for snappy_pipeline"""

import os
import sys
from collections import OrderedDict
from collections.abc import Iterator, MutableMapping
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, AnyStr, Dict

import ruamel.yaml as ruamel_yaml

from .models import SnappyModel, SnappyStepModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class SkipLibraryWarning(UserWarning):
    """Raised when libraries are skipped."""


class ConfigError(ValueError):
    """A problem of the project configuration or the sample sheets that the user has to fix.

    The message states the problem. The orchestrator attaches the task it belongs to with
    :func:`in_task` around its checks of one task, ``task_produces`` and the constructors of the
    workflow classes; ``str()`` then reads ``Task 'name' (step step): message``. The main
    Snakefile and ``snappy check-files`` print it and exit, without a traceback; in an input
    function, Snakemake names the rule and the wildcards. A ``ValueError`` or ``TypeError``
    marks a programming error.
    """

    def __init__(self, message: str):
        super().__init__(message)
        self.message = message
        self.task: TaskModel | None = None

    def __str__(self) -> str:
        if self.task is None:
            return self.message
        return f"Task {self.task.name!r} (step {self.task.step}): {self.message}"


@contextmanager
def in_task(task: "TaskModel") -> Iterator[None]:
    """Attach ``task`` to a :class:`ConfigError` raised in the block, unless it has its task."""
    try:
        yield
    except ConfigError as e:
        if e.task is None:
            e.task = task
        raise


class UnsupportedActionException(Exception):
    """Raised when user try to call action that isn't supported."""


class UnknownFiltrationSourceException(Exception):
    """Raised when user try to request an unknown filtration source."""


def validate_config[C: SnappyStepModel](
    config: dict[Any, Any],
    model: type[C],
) -> C:
    return model(**config)


def print_config(config: dict[str, Any], file=sys.stderr):
    """Print human-readable version of configuration to ``file``"""
    print("\nConfiguration", file=file)
    print("-------------\n", file=file)
    yaml = ruamel_yaml.YAML()
    return yaml.dump(config, stream=file)


if TYPE_CHECKING:
    from snappy_pipeline.workflow_model import TaskModel
    from snappy_pipeline.workflows.abstract import BaseStep


def print_sample_sheets(step: "BaseStep", file=sys.stderr):
    """Print loaded sample sheets from ``BaseStep`` in human-readable format"""
    for info in step.data_set_infos:
        print("\nSample Sheet {}".format(info.sheet_path), file=file)
        print("-------------" + "-" * len(info.sheet_path) + "\n", file=file)
        yaml = ruamel_yaml.YAML()
        return yaml.dump(info.sheet.json_data, stream=file)


type DictLike = Dict | MutableMapping | SnappyModel


def merge_dictlikes[D](dict1: DictLike, dict2: DictLike, dict_class: D = OrderedDict) -> D:
    """Merge dictionary/model ``dict2`` into ``dict1``"""

    def _merge_inner(d1: DictLike, d2: DictLike) -> D:
        DICT_LIKE = DictLike.__value__
        for k in d1.keys() | d2.keys():
            if k in d1 and k in d2:
                if isinstance(d1[k], DICT_LIKE) and isinstance(d2[k], DICT_LIKE):
                    yield k, dict_class(_merge_inner(d1[k], d2[k]))
                else:
                    # If one of the values is not a dict, you can't continue
                    # merging it.  Value from second dict overrides one in
                    # first and we move on.
                    yield k, d2[k]
            elif k in d1:
                yield k, d1[k]
            else:
                yield k, d2[k]

    return dict_class(_merge_inner(dict1, dict2))


def snakefile_path(step_name: str) -> AnyStr:
    """Return absolute path to Snakefile for the given step name"""
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "workflows", step_name, "Snakefile")
    )
