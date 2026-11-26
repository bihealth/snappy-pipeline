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

from pydantic import ValidationError
from snakemake.pathvars import Pathvars

from snappy_pipeline.base import ConfigError, in_task
from snappy_pipeline.log_archive import archive_after_run, record_git_state
from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.reads import FilesTable, read_files_tables
from snappy_pipeline.workflow_model import ConfigModel, TaskModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    ExpectedIndexedRegions,
    ExpectedPathSchema,
    ExpectedReference,
    ExpectedVcf,
    field_contract,
    select_signature,
)

if TYPE_CHECKING:
    from snakemake.api import Workflow

    from snappy_pipeline.workflows.abstract import BaseStep

logger = logging.getLogger(__name__)

#: Steps on hold, not fully ported yet: they link their result files with their own jobs and have no docs page
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
    #: Task name -> ``depends_on`` field -> the upstream signature the field reads (a tuple for a
    #: list field), picked by :func:`select_signature` from what the upstream task produces
    selected: Mapping[str, Mapping[str, DataSignature | tuple[DataSignature, ...]]] = (
        dataclass_field(default_factory=dict)
    )
    #: Task name -> the level the task writes its results at, ``library`` or ``related``
    #: (see :func:`resolve_group_by`)
    group_by: Mapping[str, str] = dataclass_field(default_factory=dict)
    #: Task name (``data_sets`` for the data sets) -> the paths in its configuration, resolved.
    file_paths: Mapping[str, tuple[str, ...]] = dataclass_field(default_factory=dict)
    #: ``include_background`` -> the library table without relationships, built once per run
    #: and shared by all tasks (``BaseStep.library_table``)
    library_tables: dict[bool, Any] = dataclass_field(default_factory=dict, compare=False)
    #: ``include_background`` -> the donor table of the germline sheets (``BaseStep.pedigree_donors``)
    donor_tables: dict[bool, Any] = dataclass_field(default_factory=dict, compare=False)
    #: Snakemake's ``pathvars`` of the configuration, expanded in the paths of files tables
    pathvars: Any = None
    #: The files tables of an ``external_data`` task (its ``library_files``) -> their contents,
    #: read once per run (:meth:`files_table`)
    _files_tables: dict[tuple[str, ...], FilesTable] = dataclass_field(
        default_factory=dict, compare=False, repr=False
    )
    #: The data sets with their sample sheets, loaded on first use (:attr:`data_set_infos`)
    _data_sets: list[Any] = dataclass_field(default_factory=list, compare=False, repr=False)

    @property
    def data_set_infos(self) -> list[Any]:
        """The data sets of the configuration with their loaded sample sheets, once per run"""
        from snappy_pipeline.library_table import DataSetInfo

        if not self._data_sets and self.model.data_sets:
            self._data_sets.extend(
                DataSetInfo(
                    name,
                    data_set.file,
                    data_set.type,
                    data_set.is_background,
                    data_set.naming_scheme,
                    data_set.pedigree_field,
                )
                for name, data_set in self.model.data_sets.items()
            )
        return self._data_sets

    def files_table(self, config: Any) -> FilesTable:
        """Return the files of the libraries that the files tables of the ``external_data`` task
        config ``config`` list, read once per run; raises :class:`ConfigError` for a bad table."""
        key = tuple(config.library_files)
        if key not in self._files_tables:
            try:
                self._files_tables[key] = read_files_tables(
                    key, config.patterns_find_reads(), config.mixed_se_pe, self.pathvars
                )
            except ValueError as e:
                raise ConfigError(str(e)) from None
        return self._files_tables[key]

    def task(self, name: str) -> TaskModel:
        """Return the task called ``name``."""
        return next(task for task in self.tasks if task.name == name)


def load_project(config: Mapping[str, Any], project_dir: str) -> Project:
    """Validate ``config`` and resolve its task graph.

    Raises :class:`ConfigError` for invalid values (the pydantic errors, with the task and step
    names), unknown steps, duplicate task names, ``depends_on`` values that name no task,
    dependency cycles, upstream tasks that do not produce what a ``depends_on`` field requires,
    ``FILTER`` IDs that could mean two things (see :func:`_check_filter_ids`), and indexes built
    from another genome than the reader's ``depends_on.reference``.
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

    try:
        model = ConfigModel.model_validate(config, context=context("data_sets"))
    except ValidationError as e:
        raise ConfigError(f"{CONFIG_FILE} is invalid:\n{_format_validation_error(e)}") from None

    names = [task.name for task in model.tasks]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ConfigError(f"Task names must be unique; duplicates: {', '.join(duplicates)}")

    task_configs: dict[str, SnappyStepModel] = {}
    dependencies: dict[str, dict[str, str]] = {}
    for task in model.tasks:
        with in_task(task):
            workflow_cls = WORKFLOW_REGISTRY.get(task.step)
            if workflow_cls is None:
                raise ConfigError(
                    f"unknown step; known steps: {', '.join(sorted(WORKFLOW_REGISTRY))}"
                )
            try:
                task_config = workflow_cls.config_model_class.model_validate(
                    task.config, context=context(task.name)
                )
            except ValidationError as e:
                raise ConfigError(
                    "invalid configuration:\n" + _format_validation_error(e)
                ) from None
            if getattr(task_config, "background_selection", None) and not (
                workflow_cls.uses_background or workflow_cls.background_results
            ):
                raise ConfigError(
                    "the step uses no background libraries, so background_selection has no "
                    "effect; remove it"
                )
            task_configs[task.name] = task_config
            dependencies[task.name] = _resolve_dependencies(task.name, task_config, names)

    tasks = tuple(_topological_order(model.tasks, dependencies))
    signatures, selected = _task_signatures(tasks, task_configs, dependencies)
    project = Project(
        config=config,
        model=model,
        project_dir=project_dir,
        config_paths=(os.path.join(project_dir, CONFIG_FILE),),
        tasks=tasks,
        task_configs=task_configs,
        dependencies=dependencies,
        signatures=signatures,
        selected=selected,
        group_by=resolve_group_by(tasks, task_configs, dependencies),
        file_paths={name: tuple(paths) for name, paths in file_paths.items()},
        pathvars=pathvars,
    )
    _check_granularity(project)
    _check_filter_ids(project)
    _check_index_references(project)
    return project


def _format_validation_error(error: ValidationError) -> str:
    """Return one line per error: the field path, the message and the offending value."""
    lines = []
    for item in error.errors(include_url=False):
        path = ".".join(str(part) for part in item["loc"])
        message = item["msg"].removeprefix("Value error, ")
        if item["type"] not in ("missing", "extra_forbidden") and not isinstance(
            item["input"], (dict, list)
        ):
            message += f" (got {item['input']!r})"
        lines.append(f"  {path}: {message}" if path else f"  {message}")
    return "\n".join(lines)


def _check_index_references(project: Project) -> None:
    """Raise if a task reads its ``depends_on.index`` from a ``reference_index`` task that builds
    it from another genome than the task's ``depends_on.reference``.

    Genomes are compared by task, and two ``external_data`` tasks that name the same FASTA file
    are one genome; a custom genome (mouse with extra contigs, human with bacteria) is one FASTA
    that both tasks name. Indexes from ``external_data`` have no genome to compare.
    """

    def genome(task_name: str) -> str:
        files = getattr(project.task_configs[task_name], "files", None) or {}
        return files.get("fasta") or task_name

    for task in project.tasks:
        with in_task(task):
            dependencies = project.dependencies[task.name]
            index, reference = dependencies.get("index"), dependencies.get("reference")
            if not reference or not index or project.task(index).step != "reference_index":
                continue
            built_from = project.dependencies[index]["reference"]
            if genome(built_from) != genome(reference):
                raise ConfigError(
                    f"the index of task {index!r} is built from the genome of task "
                    f"{built_from!r}, but depends_on.reference names task {reference!r}; "
                    "name the same genome in both"
                )


def _check_filter_ids(project: Project) -> None:
    """Raise if a ``FILTER`` value could mean two things: a task named ``PASS`` or like an ID that
    another task writes under a name of its own (variant callers, vembrane expressions, DKFZ bias
    filter; filters such as bcftools write their task name), or two tasks of one chain that write
    the same ID."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    if "PASS" in project.task_configs:
        raise ConfigError("Task name 'PASS' is reserved for the FILTER value of passing variants")
    ids = {
        name: WORKFLOW_REGISTRY[project.task(name).step].filter_ids(config)
        for name, config in project.task_configs.items()
    }
    for owner, owned in ids.items():
        clashes = sorted(owned & (set(ids) - {owner}))
        if clashes:
            raise ConfigError(
                f"Task name {clashes[0]!r} is also a FILTER ID of task {owner!r}; rename one of them"
            )
    # A VCF has the FILTER IDs of all tasks upstream of it: ID -> the task that writes it
    writers: dict[str, dict[str, str]] = {}
    for task in project.tasks:
        inherited: dict[str, str] = {}
        for upstream in upstream_tasks(project.dependencies[task.name]):
            inherited |= writers.get(upstream, {})
        shared = sorted(ids[task.name] & set(inherited))
        if shared:
            raise ConfigError(
                f"Tasks {inherited[shared[0]]!r} and {task.name!r} both write the FILTER ID "
                f"{shared[0]!r} along one chain; rename one of them"
            )
        writers[task.name] = inherited | dict.fromkeys(ids[task.name], task.name)


def resolve_group_by(
    tasks: tuple[TaskModel, ...],
    task_configs: Mapping[str, SnappyStepModel],
    dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]],
) -> dict[str, str]:
    """Return task name -> ``library`` or ``related``: the level a task writes its results at.

    A task declares its own ``group_by``, else its step's :meth:`default_group_by`, else what its
    upstream tasks declare (the first in ``depends_on`` order that declares one, transitively), so
    a filtration after joint calling writes per group of related libraries. A task that declares
    nothing writes per library, and passes nothing on, so that default does not spread downstream.
    ``tasks`` must be in topological order.
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    declared: dict[str, str | None] = {}
    for task in tasks:
        config = task_configs[task.name]
        level = getattr(config, "group_by", None)
        if level is None:
            level = WORKFLOW_REGISTRY[task.step].default_group_by(config) or None
        if level is None:
            upstreams = upstream_tasks(dependencies.get(task.name, {}))
            level = next((declared[name] for name in upstreams if declared.get(name)), None)
        declared[task.name] = level
    return {name: level or "library" for name, level in declared.items()}


def _check_granularity(project: Project) -> None:
    """Raise if a task reads results that its upstream task writes at another level.
    A ``depends_on`` field with a :class:`Granularity` names the level it needs (``ANY_LEVEL`` takes both); any other field reads the result of the task's own entity, so a task behind one that writes per related group must write per related group too."""
    for task in project.tasks:
        with in_task(task):
            config = project.task_configs[task.name]
            fields = type(config.depends_on).model_fields if project.dependencies[task.name] else {}
            for field, value in project.dependencies[task.name].items():
                required = field_contract(fields[field], config.tool_tag()).level
                for upstream in value if isinstance(value, tuple) else (value,):
                    level = project.group_by[upstream]
                    if required is not None:
                        if required.level is not None and level != required.level:
                            raise ConfigError(
                                f"depends_on.{field} needs results with group_by: {required.level}, "
                                f"but task {upstream!r} writes them with group_by: {level}"
                            )
                    elif level == "related" and project.group_by[task.name] != "related":
                        raise ConfigError(
                            f"task {upstream!r} writes one result per related group, but "
                            f"depends_on.{field} reads it per {project.group_by[task.name]}; set "
                            "group_by: related on this task"
                        )


def _exists(path: str) -> bool:
    """Whether ``path`` exists, or is the prefix of files such as a bwa index (``path.amb``, ...)."""
    return os.path.exists(path) or bool(glob.glob(glob.escape(path) + ".*"))


def missing_files(project: Project) -> dict[str, list[str]]:
    """Return task name (``data_sets`` for the data sets) -> the files of its configuration that do
    not exist: its ``ConfigPath`` values and the sidecars that the jobs reading an
    ``external_data`` task declare as inputs (the ``.fai`` and ``.dict`` of a genome, the index of
    a VCF or of an indexed BED file)."""
    paths = {name: list(found) for name, found in project.file_paths.items()}
    for name, config in project.task_configs.items():
        tables = getattr(config, "library_files", None)
        if tables and all(os.path.exists(table) for table in tables):
            paths[name] += project.files_table(config).paths()
    for consumer, dependencies in project.dependencies.items():
        config = project.task_configs[consumer]
        fields = type(config.depends_on).model_fields if dependencies else {}
        for field, upstream in dependencies.items():
            files = getattr(project.task_configs.get(upstream), "files", None) or {}
            schema = next(
                (
                    m.for_tool(config.tool_tag())
                    for m in fields[field].metadata
                    if isinstance(m, ExpectedPathSchema)
                ),
                None,
            )
            if schema is ExpectedReference and isinstance(files.get("fasta"), str):
                reference = ExpectedReference(fasta=files["fasta"])
                paths[upstream] += [reference.fai, reference.sequence_dict]
            elif schema is ExpectedVcf and isinstance(files.get("vcf"), str):
                paths[upstream].append(ExpectedVcf(vcf=files["vcf"]).index)
            elif schema is ExpectedIndexedRegions and isinstance(files.get("bed"), str):
                paths[upstream].append(ExpectedIndexedRegions(bed=files["bed"]).index)
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
    resolve through the project, so their output paths can be computed. A :class:`ConfigError`
    of a constructor gets its task (:func:`in_task`).
    """
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    _TASK_INSTANCES.clear()
    for task in project.tasks:
        if task_names is None or task.name in task_names:
            with in_task(task):
                _TASK_INSTANCES[task.name] = WORKFLOW_REGISTRY[task.step](
                    workflow, project, task.name
                )
    return dict(_TASK_INSTANCES)


def task_instance(task_name: str) -> BaseStep:
    """Return the workflow object of ``task_name``; called from the step Snakefiles."""
    return _TASK_INSTANCES[task_name]


def register_hooks(workflow: Workflow) -> None:
    """Record the git state of the configuration when the run starts, archive the logs of the run
    and print a banner when the workflow fails or succeeds."""

    def finish(snakemake_logs: list[str], message: str) -> None:
        archive_after_run(Path.cwd(), snakemake_logs)
        line = "*" * len(message)
        print(f"\n{line}\n{message}\n{line}\n", file=sys.stderr)

    workflow.onstart(lambda logs: record_git_state(Path.cwd(), logs))
    workflow.onerror(lambda logs: finish(logs, "Oh no! Something went wrong."))
    workflow.onsuccess(lambda logs: finish(logs, "All done; have a nice day!"))


def select_target_tasks(project: Project, target_task: str | None = None) -> list[str]:
    """Return the names of the tasks that ``snappy run`` builds and links the result files of, in config order.

    Every task, or only ``target_task`` if given. Raises :class:`ConfigError` if no task has that name.
    """
    names = [task.name for task in project.model.tasks]
    if target_task is not None:
        if target_task not in names:
            raise ConfigError(f"Unknown task {target_task!r}; configured tasks: {', '.join(names)}")
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

    """
    depends_on = getattr(task_config, "depends_on", None)
    if depends_on is None:
        return {}
    result: dict[str, str | tuple[str, ...]] = {}
    for field, value in depends_on.model_dump().items():
        if not value:
            continue
        upstreams = value if isinstance(value, list) else [value]
        for upstream in upstreams:
            if upstream == task_name:
                raise ConfigError(f"depends_on.{field} names the task itself")
            if upstream not in names:
                verb = "contains" if isinstance(value, list) else "is"
                raise ConfigError(
                    f"depends_on.{field} {verb} {upstream!r}, which is not a task; "
                    f"configured tasks: {', '.join(names)}"
                )
        if len(set(upstreams)) != len(upstreams):
            raise ConfigError(f"depends_on.{field} names a task twice")
        result[field] = tuple(upstreams) if isinstance(value, list) else value
    return result


Selected = dict[str, DataSignature | tuple[DataSignature, ...]]


def _task_signatures(
    tasks: tuple[TaskModel, ...],
    task_configs: Mapping[str, SnappyStepModel],
    dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]],
) -> tuple[dict[str, tuple[DataSignature, ...]], dict[str, Selected]]:
    """Return task name -> produced signatures and task name -> the signatures its ``depends_on``
    fields read, checking every ``depends_on`` requirement.

    Each task's ``task_produces`` receives, per field, the upstream signature it reads
    (``select_signature``); for a list field, a tuple of them. The steps read the same selection
    through ``BaseStep.get_upstream_paths``.

    ``tasks`` must be in topological order, so the upstream signatures exist when a task needs
    them.
    """
    signatures: dict[str, tuple[DataSignature, ...]] = {}
    selected: dict[str, Selected] = {}
    steps = {task.name: task.step for task in tasks}
    for task in tasks:
        with in_task(task):
            signatures[task.name], selected[task.name] = _task_signature(
                task, task_configs, dependencies, signatures, steps
            )
    return signatures, selected


def _task_signature(
    task: TaskModel,
    task_configs: Mapping[str, SnappyStepModel],
    dependencies: Mapping[str, Mapping[str, str | tuple[str, ...]]],
    signatures: Mapping[str, tuple[DataSignature, ...]],
    steps: Mapping[str, str],
) -> tuple[tuple[DataSignature, ...], Selected]:
    """Return the signatures that ``task`` produces and those its fields read, given those of
    its upstream tasks."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    config = task_configs[task.name]
    fields = type(config.depends_on).model_fields if dependencies[task.name] else {}
    upstream_signatures: Selected = {}
    for field, value in dependencies[task.name].items():
        contract = field_contract(fields[field], config.tool_tag())
        selected_all = []
        for upstream in value if isinstance(value, tuple) else (value,):
            produced = signatures[upstream]
            selected = select_signature(produced, contract.signature)
            if selected is None:
                raise ConfigError(
                    f"depends_on.{field} requires {contract.signature}, but task {upstream!r} "
                    f"produces {', '.join(map(str, produced)) or 'nothing'}"
                )
            selected_all.append(selected)
            _check_output_keys(field, contract.schema, upstream, selected, steps, task_configs)
        upstream_signatures[field] = (
            tuple(selected_all) if isinstance(value, tuple) else selected_all[0]
        )
    produces = WORKFLOW_REGISTRY[task.step].task_produces(config, upstream_signatures)
    return produces, upstream_signatures


def _check_output_keys(field, schema, upstream, signature, steps, task_configs) -> None:
    """Raise if the upstream task does not provide the keys of the field's ``ExpectedPathSchema``
    for the signature the field selected."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    if schema is None:
        return
    keys = WORKFLOW_REGISTRY[steps[upstream]].output_keys(task_configs[upstream], signature)
    if keys is None:  # a files table that does not exist yet; checked at setup
        return
    missing = [
        name for name, f in schema.model_fields.items() if f.is_required() and name not in keys
    ]
    if missing:
        raise ConfigError(
            f"depends_on.{field} needs {', '.join(missing)}, which task {upstream!r} does not "
            f"provide; it provides {', '.join(sorted(keys))}"
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
            raise ConfigError(f"Dependency cycle between tasks: {' -> '.join(cycle)}")
        state[name] = "visiting"
        for upstream in upstream_tasks(dependencies[name]):
            visit(upstream, chain + [name])
        state[name] = "done"
        ordered.append(by_name[name])

    for task in tasks:
        visit(task.name, [])
    return ordered
