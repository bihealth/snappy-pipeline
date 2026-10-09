"""Loading a project and its task graph for the orchestrator Snakefile (``snappy_pipeline/Snakefile``)."""

from __future__ import annotations

import logging
import os
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.workflow_model import ConfigModel, TaskModel
from snappy_pipeline.workflows.abstract.protocol import DATA_SETS, DataSignature, select_signature

if TYPE_CHECKING:
    from snakemake.api import Workflow

    from snappy_pipeline.workflows.abstract import BaseStep

logger = logging.getLogger(__name__)

#: Task name -> workflow object, filled by the orchestrator before it loads the step modules.
_TASK_INSTANCES: dict[str, BaseStep] = {}

#: Name of the project configuration file in the project directory.
CONFIG_FILE = "config.yaml"


@dataclass(frozen=True)
class Project:
    """A validated project configuration and its resolved task graph."""

    #: The configuration as read from ``config.yaml``, without ``null`` values.
    config: Mapping[str, Any]
    #: The validated global configuration.
    model: ConfigModel
    #: The project directory.
    work_dir: str
    #: Directories that relative paths in the configuration are resolved against.
    lookup_paths: tuple[str, ...]
    #: Configuration files the project was loaded from.
    config_paths: tuple[str, ...]
    #: The tasks, each after the tasks it depends on.
    tasks: tuple[TaskModel, ...]
    #: Task name -> validated step configuration.
    task_configs: Mapping[str, SnappyStepModel]
    #: Task name -> ``depends_on`` field -> upstream task name, for set fields only.
    dependencies: Mapping[str, Mapping[str, str]]
    #: Task name -> the DataSignatures the task produces.
    signatures: Mapping[str, tuple[DataSignature, ...]]

    def task(self, name: str) -> TaskModel:
        """Return the task called ``name``."""
        return next(task for task in self.tasks if task.name == name)


def load_project(config: Mapping[str, Any], work_dir: str) -> Project:
    """Validate ``config`` and resolve its task graph.

    Raises ``ValueError`` (or a pydantic ``ValidationError``) for unknown steps, duplicate task
    names, ``depends_on`` values that name no task, dependency cycles, and upstream tasks that
    do not produce what a ``depends_on`` field requires.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    config = _without_nulls(config)
    model = ConfigModel(**config)
    work_dir = os.path.abspath(work_dir)
    lookup_paths = (work_dir, os.path.dirname(work_dir))

    names = [task.name for task in model.tasks]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Task names must be unique; duplicates: {', '.join(duplicates)}")
    if DATA_SETS in names:
        raise ValueError(f"Task name {DATA_SETS!r} is reserved for depends_on.reads")

    task_configs: dict[str, SnappyStepModel] = {}
    dependencies: dict[str, dict[str, str]] = {}
    for task in model.tasks:
        workflow_cls = WORKFLOW_REGISTRY.get(task.step)
        if workflow_cls is None:
            raise ValueError(
                f"Task {task.name!r} uses unknown step {task.step!r}; "
                f"known steps: {', '.join(sorted(WORKFLOW_REGISTRY))}"
            )
        task_config = workflow_cls.config_model_class.model_validate(
            task.config, context={"config_lookup_paths": lookup_paths}
        )
        task_configs[task.name] = task_config
        dependencies[task.name] = _resolve_dependencies(task.name, task_config, names)

    tasks = tuple(_topological_order(model.tasks, dependencies))
    return Project(
        config=config,
        model=model,
        work_dir=work_dir,
        lookup_paths=lookup_paths,
        config_paths=(os.path.join(work_dir, CONFIG_FILE),),
        tasks=tasks,
        task_configs=task_configs,
        dependencies=dependencies,
        signatures=_task_signatures(tasks, task_configs, dependencies),
    )


def create_task_instances(workflow: Workflow, project: Project) -> dict[str, BaseStep]:
    """Create the workflow object of every task once and register it for the step Snakefiles."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    _TASK_INSTANCES.clear()
    for task in project.tasks:
        _TASK_INSTANCES[task.name] = WORKFLOW_REGISTRY[task.step](workflow, project, task.name)
    return dict(_TASK_INSTANCES)


def task_instance(task_name: str) -> BaseStep:
    """Return the workflow object of ``task_name``; called from the step Snakefiles."""
    return _TASK_INSTANCES[task_name]


def register_hooks(workflow: Workflow) -> None:
    """Print a banner when the workflow fails or succeeds."""

    def banner(message: str) -> None:
        line = "*" * len(message)
        print(f"\n{line}\n{message}\n{line}\n", file=sys.stderr)

    workflow.onerror(lambda _: banner("Oh no! Something went wrong."))
    workflow.onsuccess(lambda _: banner("All done; have a nice day!"))


def select_target_tasks(
    project: Project, target_task: str | None = None, all_tasks: bool = False
) -> list[str]:
    """Return the names of the tasks whose result files ``snappy run`` targets, in config order.

    * ``target_task`` given: only that task. Raises ``ValueError`` if no task has that name.
    * ``all_tasks``: every task.
    * Otherwise the leaf tasks: tasks that no other task depends on.
    """
    names = [task.name for task in project.model.tasks]
    if target_task is not None:
        if target_task not in names:
            raise ValueError(f"Unknown task {target_task!r}; configured tasks: {', '.join(names)}")
        return [target_task]
    if all_tasks:
        return names
    depended_on = {name for deps in project.dependencies.values() for name in deps.values()}
    return [name for name in names if name not in depended_on]


def _without_nulls(data: Any, path: str = "") -> Any:
    """Return ``data`` without ``null`` values, so that the model defaults apply to them."""
    if isinstance(data, Mapping):
        result = {}
        for key, value in data.items():
            key_path = f"{path}.{key}" if path else str(key)
            if value is None:
                logger.info("Configuration key %s has no value; using its default", key_path)
            else:
                result[key] = _without_nulls(value, key_path)
        return result
    if isinstance(data, list):
        return [_without_nulls(value, f"{path}[{i}]") for i, value in enumerate(data)]
    return data


def _resolve_dependencies(
    task_name: str, task_config: SnappyStepModel, names: list[str]
) -> dict[str, str]:
    """Return ``depends_on`` field -> upstream task name for the set fields of one task.

    ``reads: data_sets`` names no task, so it is not a dependency.
    """
    depends_on = getattr(task_config, "depends_on", None)
    if depends_on is None:
        return {}
    result = {}
    for field, upstream in depends_on.model_dump().items():
        if not upstream or (field == "reads" and upstream == DATA_SETS):
            continue
        if upstream == task_name:
            raise ValueError(f"Task {task_name!r}: depends_on.{field} names the task itself")
        if upstream not in names:
            raise ValueError(
                f"Task {task_name!r}: depends_on.{field} is {upstream!r}, which is not a task; "
                f"configured tasks: {', '.join(names)}"
            )
        result[field] = upstream
    return result


def _task_signatures(
    tasks: tuple[TaskModel, ...],
    task_configs: Mapping[str, SnappyStepModel],
    dependencies: Mapping[str, Mapping[str, str]],
) -> dict[str, tuple[DataSignature, ...]]:
    """Return task name -> produced signatures, checking every ``depends_on`` requirement.

    Each task's ``task_produces`` receives, per field, the upstream signature it reads
    (``select_signature``).

    ``tasks`` must be in topological order, so the upstream signatures exist when a task needs
    them.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    signatures: dict[str, tuple[DataSignature, ...]] = {}
    for task in tasks:
        config = task_configs[task.name]
        fields = type(config.depends_on).model_fields if dependencies[task.name] else {}
        upstream_signatures = {}
        for field, upstream in dependencies[task.name].items():
            required = next(
                (m for m in fields[field].metadata if isinstance(m, DataSignature)), None
            )
            produced = signatures[upstream]
            selected = select_signature(produced, required)
            if selected is None:
                raise ValueError(
                    f"Task {task.name!r}: depends_on.{field} requires {required}, but task "
                    f"{upstream!r} produces {', '.join(map(str, produced)) or 'nothing'}"
                )
            upstream_signatures[field] = selected
        signatures[task.name] = WORKFLOW_REGISTRY[task.step].task_produces(
            config, upstream_signatures
        )
    return signatures


def _topological_order(
    tasks: list[TaskModel], dependencies: Mapping[str, Mapping[str, str]]
) -> list[TaskModel]:
    """Return ``tasks`` with every task after its dependencies, keeping config order otherwise."""
    by_name = {task.name: task for task in tasks}
    ordered: list[TaskModel] = []
    state: dict[str, str] = {}  # "visiting" or "done"

    def visit(name: str, chain: list[str]) -> None:
        if state.get(name) == "done":
            return
        if state.get(name) == "visiting":
            cycle = chain[chain.index(name) :] + [name]
            raise ValueError(f"Dependency cycle between tasks: {' -> '.join(cycle)}")
        state[name] = "visiting"
        for upstream in dependencies[name].values():
            visit(upstream, chain + [name])
        state[name] = "done"
        ordered.append(by_name[name])

    for task in tasks:
        visit(task.name, [])
    return ordered
