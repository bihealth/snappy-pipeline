"""The rule names of a project: ``<task>_<part>_<action>`` for the jobs of a step part and
``<task>_results`` for the job that links a task's results."""

from collections.abc import Collection, Mapping
from typing import Any, NamedTuple

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


def parse_rule_name(name: str, parts: TaskParts) -> RuleName | None:
    """The task, part and action of a rule, None for a rule of no task (``all``).

    A task name may be the prefix of another one's (``mapping``, ``mapping_qc``), so the task
    whose parts match the rest of the name wins (``mapping_qc_summary_run`` is the ``qc_summary``
    part of ``mapping``); without a match, the longest task name does.
    """
    candidates = sorted((t for t in parts if name.startswith(f"{t}_")), key=len, reverse=True)
    for task in candidates:
        rest = name.removeprefix(f"{task}_")
        for part, actions in sorted(parts[task].items(), key=lambda item: -len(item[0])):
            action = rest.removeprefix(f"{part}_")
            if action != rest and action in actions:
                return RuleName(task, part, action)
    if candidates:
        return RuleName(candidates[0], "", "")
    return None
