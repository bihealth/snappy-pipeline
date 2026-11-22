# -*- coding: utf-8 -*-
"""What the task questions of ``snappy init`` and ``snappy task add`` read from the step models:
the tools of a step, the ``depends_on`` keys a tool requires, which of the tasks added so far
can provide a key, which steps could, and the settings of an ``external_data`` task for a key.
"""

from dataclasses import dataclass
from typing import Any

from snappy_pipeline.models import RequiredBy, _is_optional, reader_markers, tool_values
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedAlignments,
    ExpectedKits,
    ExpectedVariants,
    FieldContract,
    field_contract,
)

#: The step of existing data, which can provide any key
EXTERNAL = "external_data"
#: The FASTQ source that is no task
DATA_SETS = "data_sets"
#: Schemas of files per library, which an ``external_data`` task gives through search patterns
PER_LIBRARY = (ExpectedAlignments, ExpectedVariants)


@dataclass(frozen=True)
class TaskEntry:
    """A task of the project: name, step, tool and, for an ``external_data`` task, its data."""

    name: str
    step: str
    tool: str | None = None
    produces: DataSignature | None = None


@dataclass(frozen=True)
class Key:
    """A ``depends_on`` key a task must set, with its contract and the field's first sentence."""

    name: str
    contract: FieldContract
    description: str


def entries_of(config_yaml: dict | None) -> list[TaskEntry]:
    """The tasks of ``config.yaml`` in order; a task without name or step is left out."""
    entries = []
    for item in (config_yaml or {}).get("tasks") or []:
        if not isinstance(item, dict) or not item.get("name") or not item.get("step"):
            continue
        config = item.get("config") or {}
        produces = None
        declared = config.get("produces")
        if item["step"] == EXTERNAL and isinstance(declared, dict) and declared.get("type"):
            produces = DataSignature(
                DataType(declared["type"]), frozenset(declared.get("tags") or [])
            )
        entries.append(
            TaskEntry(str(item["name"]), str(item["step"]), config.get("tool"), produces)
        )
    return entries


def step_tools(step: str) -> list[str]:
    """The tool values of ``step``; empty for a step without a tool."""
    return tool_values(WORKFLOW_REGISTRY[step].config_model_class)


def required_keys(step: str, tool: str | None) -> list[Key]:
    """The ``depends_on`` keys a task of ``step`` with ``tool`` must set: those without a default
    and those a ``RequiredBy`` marker of the tool needs regardless of a setting."""
    field = WORKFLOW_REGISTRY[step].config_model_class.model_fields.get("depends_on")
    if field is None:
        return []
    keys = []
    for name, sub_field in field.annotation.model_fields.items():
        needed = [
            marker
            for marker in reader_markers(sub_field)
            if isinstance(marker, RequiredBy)
            and marker.when is None
            and (not marker.tools or tool in marker.tools)
        ]
        if not _is_optional(sub_field) or needed:
            description = (sub_field.description or "").split("\n")[0]
            keys.append(Key(name, field_contract(sub_field, tool), description))
    return keys


def _produces_type(step: str, data_type: DataType) -> bool:
    return any(signature.type == data_type for signature in WORKFLOW_REGISTRY[step].produces)


def candidates(key: Key, entries: list[TaskEntry]) -> list[str]:
    """The names of the tasks in ``entries`` that can provide ``key``, newest first:
    an ``external_data`` task by what it produces, another task by the type its step produces;
    ``data_sets`` first for the reads."""
    names = [DATA_SETS] if key.name == "reads" else []
    signature = key.contract.signature
    for entry in reversed(entries):
        if signature is None:
            fits = True
        elif entry.step == EXTERNAL:
            fits = entry.produces is not None and entry.produces.satisfies(signature)
        else:
            fits = _produces_type(entry.step, signature.type)
        if fits:
            names.append(entry.name)
    return names


def providers(key: Key) -> list[str]:
    """The steps a task can have to provide ``key``: ``external_data``, then the steps that
    produce the type of its requirement."""
    signature = key.contract.signature
    steps = sorted(
        step
        for step in WORKFLOW_REGISTRY
        if step != EXTERNAL and (signature is None or _produces_type(step, signature.type))
    )
    return [EXTERNAL, *steps]


def external_data_settings(key: Key) -> dict[str, Any]:
    """The settings of an ``external_data`` task that provides ``key``: ``produces`` with the type
    and the plain tags of the requirement, and the expected file keys without paths, to fill in
    (``kits`` as a whole for kit files; nothing for files per library, which search patterns
    find)."""
    settings: dict[str, Any] = {}
    signature = key.contract.signature
    if signature is not None:
        tags = sorted(
            tag[0] if isinstance(tag, tuple) else tag
            for tag in signature.tags
            if isinstance(tag, tuple) or not tag.startswith("-")
        )
        settings["produces"] = {"type": signature.type.value, "tags": tags}
    schema = key.contract.schema
    if schema is ExpectedKits:
        settings["kits"] = None
    elif schema is not None and schema not in PER_LIBRARY:
        settings["files"] = {name: None for name in schema.model_fields}
    return settings


def entry_for(name: str, step: str, tool: str | None, settings: dict[str, Any]) -> TaskEntry:
    """The entry of a task written with ``settings`` (see :func:`external_data_settings`)."""
    produces = settings.get("produces")
    signature = None
    if produces:
        signature = DataSignature(DataType(produces["type"]), frozenset(produces["tags"]))
    return TaskEntry(name, step, tool, signature)
