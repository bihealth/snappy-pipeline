"""Helpers for the orchestrator Snakefile (``snappy_pipeline/Snakefile``)."""

from collections.abc import Mapping, Sequence
from typing import Any


def select_target_tasks(
    tasks: Sequence[Mapping[str, Any]], target_task: str | None = None, all_tasks: bool = False
) -> list[str]:
    """Return the names of the tasks whose result files ``snappy run`` targets, in config order.

    * ``target_task`` given: only that task. Raises ``ValueError`` if no task has that name.
    * ``all_tasks``: every task.
    * Otherwise the leaf tasks: tasks that no other task names in its ``config.depends_on``.
      Empty values (unset optional dependencies) do not count.
    """
    names = [task["name"] for task in tasks]
    if target_task is not None:
        if target_task not in names:
            raise ValueError(f"Unknown task {target_task!r}; configured tasks: {', '.join(names)}")
        return [target_task]
    if all_tasks:
        return names
    depended_on = {
        value
        for task in tasks
        for value in ((task.get("config") or {}).get("depends_on") or {}).values()
        if value
    }
    return [name for name in names if name not in depended_on]
