# -*- coding: utf-8 -*-
"""Basic utility code for snappy_pipeline"""

from collections.abc import Iterator
from contextlib import contextmanager
from typing import TYPE_CHECKING

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

if TYPE_CHECKING:
    from snappy_pipeline.workflow_model import TaskModel

#: Name of the project configuration file in the project directory.
CONFIG_FILE = "config.yaml"


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
