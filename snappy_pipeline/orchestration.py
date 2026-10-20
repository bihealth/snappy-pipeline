"""Loading a project and its task graph for the orchestrator Snakefile (``snappy_pipeline/Snakefile``)."""

from __future__ import annotations

import glob
import logging
import os
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from snakemake.pathvars import Pathvars

from snappy_pipeline.log_archive import archive_after_run
from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.reads import ReadDiscovery
from snappy_pipeline.workflow_model import ConfigModel, TaskModel
from snappy_pipeline.workflows.abstract.protocol import (
    DATA_SETS,
    DataSignature,
    ExpectedPathSchema,
    ExpectedReference,
    select_signature,
)

if TYPE_CHECKING:
    from snakemake.api import Workflow

    from snappy_pipeline.workflows.abstract import BaseStep

logger = logging.getLogger(__name__)

#: Steps on hold, not fully ported yet: they publish with their own link jobs and have no docs page
ON_HOLD_STEPS = frozenset({"create_proteome"})

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
    #: The project directory: the directory of ``config.yaml``, which relative paths in the
    #: configuration are resolved against.
    project_dir: str
    #: Configuration files the project was loaded from.
    config_paths: tuple[str, ...]
    #: The tasks, each after the tasks it depends on.
    tasks: tuple[TaskModel, ...]
    #: Task name -> validated step configuration.
    task_configs: Mapping[str, SnappyStepModel]
    #: Task name -> ``depends_on`` field -> upstream task name, or a tuple of them for a field
    #: that takes a list, for set fields only.
    dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]]
    #: Task name -> the DataSignatures the task produces.
    signatures: Mapping[str, tuple[DataSignature, ...]]
    #: Task name (``data_sets`` for the data sets) -> the paths in its configuration, resolved.
    file_paths: Mapping[str, tuple[str, ...]] = dataclass_field(default_factory=dict)
    #: FASTQ discovery shared by all tasks, so each search path is walked once per run
    read_discovery: ReadDiscovery = dataclass_field(default_factory=ReadDiscovery, compare=False)
    #: ``include_background`` -> the library table without relationships, built once per run
    #: and shared by all tasks (``BaseStep.build_library_dataframe``)
    library_tables: dict[bool, Any] = dataclass_field(default_factory=dict, compare=False)

    def task(self, name: str) -> TaskModel:
        """Return the task called ``name``."""
        return next(task for task in self.tasks if task.name == name)


def load_project(config: Mapping[str, Any], project_dir: str) -> Project:
    """Validate ``config`` and resolve its task graph.

    Raises ``ValueError`` (or a pydantic ``ValidationError``) for unknown steps, duplicate task
    names, ``depends_on`` values that name no task, dependency cycles, and upstream tasks that
    do not produce what a ``depends_on`` field requires.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    config = _without_nulls(config)
    project_dir = os.path.abspath(project_dir)
    # Snakemake's defaults and semantics, so that a path means the same in snappy and in rules
    pathvars = Pathvars.with_defaults()
    pathvars.update(Pathvars.from_config(config))
    file_paths: dict[str, tuple[str, ...]] = {}

    def context(name: str) -> dict[str, Any]:
        paths: list[str] = []
        file_paths[name] = paths  # filled while validating
        return {"project_dir": project_dir, "pathvars": pathvars, "paths": paths}

    model = ConfigModel.model_validate(config, context=context(DATA_SETS))

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
            task.config, context=context(task.name)
        )
        task_configs[task.name] = task_config
        dependencies[task.name] = _resolve_dependencies(task.name, task_config, names)

    tasks = tuple(_topological_order(model.tasks, dependencies))
    return Project(
        config=config,
        model=model,
        project_dir=project_dir,
        config_paths=(os.path.join(project_dir, CONFIG_FILE),),
        tasks=tasks,
        task_configs=task_configs,
        dependencies=dependencies,
        signatures=_task_signatures(tasks, task_configs, dependencies),
        file_paths={name: tuple(paths) for name, paths in file_paths.items()},
    )


def _exists(path: str) -> bool:
    """Whether ``path`` exists, or is the prefix of files such as a bwa index (``path.amb``, ...)."""
    return os.path.exists(path) or bool(glob.glob(glob.escape(path) + ".*"))


def missing_files(project: Project) -> dict[str, list[str]]:
    """Return task name (``data_sets`` for the data sets) -> the files of its configuration that do
    not exist: its ``ConfigPath`` values and, for an ``external_data`` genome, the ``.fai`` and
    ``.dict`` that the jobs reading it declare as inputs."""
    paths = {name: list(found) for name, found in project.file_paths.items()}
    for dependencies in project.dependencies.values():
        genome = dependencies.get("reference")
        files = getattr(project.task_configs.get(genome), "files", None) or {}
        if isinstance(files.get("fasta"), str):
            reference = ExpectedReference(fasta=files["fasta"])
            paths[genome] += [reference.fai, reference.sequence_dict]
    missing = {
        name: [path for path in dict.fromkeys(found) if not _exists(path)]
        for name, found in paths.items()
    }
    return {name: found for name, found in missing.items() if found}


def create_task_instances(
    workflow: Workflow, project: Project, task_names: list[str] | None = None
) -> dict[str, BaseStep]:
    """Create the workflow object of every task once and register it for the step Snakefiles.

    ``task_names`` limits this to some tasks (``snappy run --frozen``); upstream tasks still
    resolve through the project, so their output paths can be computed.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    _TASK_INSTANCES.clear()
    for task in project.tasks:
        if task_names is None or task.name in task_names:
            _TASK_INSTANCES[task.name] = WORKFLOW_REGISTRY[task.step](workflow, project, task.name)
    return dict(_TASK_INSTANCES)


def task_instance(task_name: str) -> BaseStep:
    """Return the workflow object of ``task_name``; called from the step Snakefiles."""
    return _TASK_INSTANCES[task_name]


def register_hooks(workflow: Workflow) -> None:
    """Archive the logs of the run and print a banner when the workflow fails or succeeds."""

    def finish(snakemake_logs: list[str], message: str) -> None:
        archive_after_run(Path.cwd(), snakemake_logs)
        line = "*" * len(message)
        print(f"\n{line}\n{message}\n{line}\n", file=sys.stderr)

    workflow.onerror(lambda logs: finish(logs, "Oh no! Something went wrong."))
    workflow.onsuccess(lambda logs: finish(logs, "All done; have a nice day!"))


def select_target_tasks(project: Project, target_task: str | None = None) -> list[str]:
    """Return the names of the tasks that ``snappy run`` builds and publishes, in config order.

    Every task, or only ``target_task`` if given. Raises ``ValueError`` if no task has that name.
    """
    names = [task.name for task in project.model.tasks]
    if target_task is not None:
        if target_task not in names:
            raise ValueError(f"Unknown task {target_task!r}; configured tasks: {', '.join(names)}")
        return [target_task]
    return names


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


def upstream_tasks(dependencies: Mapping[str, str | tuple[str, ...]]) -> list[str]:
    """Return the upstream task names of one task's ``dependencies``, in field order."""
    return [
        name
        for upstream in dependencies.values()
        for name in (upstream if isinstance(upstream, tuple) else (upstream,))
    ]


def _resolve_dependencies(
    task_name: str, task_config: SnappyStepModel, names: list[str]
) -> dict[str, str | tuple[str, ...]]:
    """Return ``depends_on`` field -> upstream task name (a tuple for a list field) for the set
    fields of one task.

    ``reads: data_sets`` names no task, so it is not a dependency.
    """
    depends_on = getattr(task_config, "depends_on", None)
    if depends_on is None:
        return {}
    result: dict[str, str | tuple[str, ...]] = {}
    for field, value in depends_on.model_dump().items():
        if not value or (field == "reads" and value == DATA_SETS):
            continue
        upstreams = value if isinstance(value, list) else [value]
        for upstream in upstreams:
            if upstream == task_name:
                raise ValueError(f"Task {task_name!r}: depends_on.{field} names the task itself")
            if upstream not in names:
                verb = "contains" if isinstance(value, list) else "is"
                raise ValueError(
                    f"Task {task_name!r}: depends_on.{field} {verb} {upstream!r}, which is not a "
                    f"task; configured tasks: {', '.join(names)}"
                )
        if len(set(upstreams)) != len(upstreams):
            raise ValueError(f"Task {task_name!r}: depends_on.{field} names a task twice")
        result[field] = tuple(upstreams) if isinstance(value, list) else value
    return result


def _task_signatures(
    tasks: tuple[TaskModel, ...],
    task_configs: Mapping[str, SnappyStepModel],
    dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]],
) -> dict[str, tuple[DataSignature, ...]]:
    """Return task name -> produced signatures, checking every ``depends_on`` requirement.

    Each task's ``task_produces`` receives, per field, the upstream signature it reads
    (``select_signature``); for a list field, a tuple of them.

    ``tasks`` must be in topological order, so the upstream signatures exist when a task needs
    them.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    signatures: dict[str, tuple[DataSignature, ...]] = {}
    steps = {task.name: task.step for task in tasks}
    for task in tasks:
        config = task_configs[task.name]
        fields = type(config.depends_on).model_fields if dependencies[task.name] else {}
        upstream_signatures = {}
        for field, value in dependencies[task.name].items():
            required = next(
                (m for m in fields[field].metadata if isinstance(m, DataSignature)), None
            )
            if required is not None:
                required = required.for_tool(config.tool_tag())
            selected_all = []
            for upstream in value if isinstance(value, tuple) else (value,):
                produced = signatures[upstream]
                selected = select_signature(produced, required)
                if selected is None:
                    raise ValueError(
                        f"Task {task.name!r}: depends_on.{field} requires {required}, but task "
                        f"{upstream!r} produces {', '.join(map(str, produced)) or 'nothing'}"
                    )
                selected_all.append(selected)
                _check_output_keys(
                    task.name, field, fields[field], upstream, selected, steps, task_configs
                )
            upstream_signatures[field] = (
                tuple(selected_all) if isinstance(value, tuple) else selected_all[0]
            )
        signatures[task.name] = WORKFLOW_REGISTRY[task.step].task_produces(
            config, upstream_signatures
        )
    return signatures


def _check_output_keys(
    task_name, field, field_info, upstream, signature, steps, task_configs
) -> None:
    """Raise if the upstream task does not provide the keys of the field's ``ExpectedPathSchema``
    for the signature the field selected."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    tool = task_configs[task_name].tool_tag()
    schema = next(
        (m.for_tool(tool) for m in field_info.metadata if isinstance(m, ExpectedPathSchema)),
        None,
    )
    if schema is None:
        return
    keys = WORKFLOW_REGISTRY[steps[upstream]].output_keys(task_configs[upstream], signature)
    missing = [
        name for name, f in schema.model_fields.items() if f.is_required() and name not in keys
    ]
    if missing:
        raise ValueError(
            f"Task {task_name!r}: depends_on.{field} needs {', '.join(missing)}, which task "
            f"{upstream!r} does not provide; it provides {', '.join(sorted(keys))}"
        )


def _topological_order(
    tasks: list[TaskModel], dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]]
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
        for upstream in upstream_tasks(dependencies[name]):
            visit(upstream, chain + [name])
        state[name] = "done"
        ordered.append(by_name[name])

    for task in tasks:
        visit(task.name, [])
    return ordered
