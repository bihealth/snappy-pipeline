# -*- coding: utf-8 -*-
"""Tests for target selection in ``snappy_pipeline.orchestration``."""

import pytest

from snappy_pipeline.orchestration import select_target_tasks


def _task(name, **depends_on):
    return {"step": "dummy", "name": name, "config": {"depends_on": depends_on}}


#: mapping -> calling -> annotation -> filtration, plus a QC task on the mapping.
TASKS = [
    _task("mapping"),
    _task("calling", ngs_mapping="mapping"),
    _task("annotation", variant="calling"),
    _task("filtration", variant="annotation", ngs_mapping="mapping"),
    _task("qc", ngs_mapping="mapping"),
]


def test_default_targets_leaf_tasks_in_config_order():
    assert select_target_tasks(TASKS) == ["filtration", "qc"]


def test_unset_optional_dependencies_do_not_count():
    tasks = [_task("mapping"), _task("calling", ngs_mapping="mapping", adapter_trimming="")]
    assert select_target_tasks(tasks) == ["calling"]


def test_tasks_without_depends_on_are_leaves_unless_depended_on():
    tasks = [{"step": "dummy", "name": "a"}, {"step": "dummy", "name": "b", "config": {}}]
    assert select_target_tasks(tasks) == ["a", "b"]


def test_single_task():
    assert select_target_tasks(TASKS, target_task="calling") == ["calling"]


def test_unknown_single_task_raises():
    with pytest.raises(ValueError, match="Unknown task 'calls'"):
        select_target_tasks(TASKS, target_task="calls")


def test_all_tasks():
    assert select_target_tasks(TASKS, all_tasks=True) == [t["name"] for t in TASKS]
