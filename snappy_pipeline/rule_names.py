"""The rule names of a project: ``<task>__<part>_<action>`` for the jobs of a step part and
``<task>__results`` for the job that links a task's results."""

from collections.abc import Collection, Mapping
from typing import Any, NamedTuple

#: Between the task name and the name of the rule in its step; a task name must not contain it
TASK_SEPARATOR = "__"

#: Task name -> part name -> the part's actions
TaskParts = Mapping[str, Mapping[str, Collection[str]]]


class RuleName(NamedTuple):
    task: str
    part: str  # empty for the results job or a rule outside the pattern
    action: str


def task_parts(instances: Mapping[str, Any]) -> dict[str, dict[str, tuple[str, ...]]]:
    """The ``TaskParts`` of the step instances by task name."""
    return {
        task: {part.name: tuple(part.actions) for part in step.sub_steps.values()}
        for task, step in instances.items()
    }


def split_task(name: str) -> tuple[str, str] | None:
    """The task and the rest of a rule name, None for a rule of no task (``all``)."""
    task, separator, rest = name.partition(TASK_SEPARATOR)
    return (task, rest) if separator else None


def parse_rule_name(name: str, parts: TaskParts) -> RuleName | None:
    """The task, part and action of a rule, None for a rule of no task (``all``).

    The part and action stay empty for the results job and a rule outside the pattern.
    """
    split = split_task(name)
    if split is None:
        return None
    task, rest = split
    for part, actions in sorted(parts.get(task, {}).items(), key=lambda item: -len(item[0])):
        action = rest.removeprefix(f"{part}_")
        if action != rest and action in actions:
            return RuleName(task, part, action)
    return RuleName(task, "", "")
