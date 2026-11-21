# -*- coding: utf-8 -*-
"""Base classes for the actual pipeline steps"""

import logging
import os
import os.path
import inspect
import math
import re
import tempfile
import typing
from collections.abc import Mapping
from dataclasses import dataclass, field
from fnmatch import fnmatch
from functools import cached_property, wraps
from typing import Any, Callable

import pydantic
from snakemake.api import Workflow
from snakemake.io import expand, report, temp
from snakemake.iocontainers import InputFiles, OutputFiles, Wildcards

from snappy_pipeline.base import (
    ConfigError,
    UnsupportedActionException,
)
from snappy_pipeline.models import (
    Intermediate,
    LibraryFolders,
    RelationshipDefinition,
    SnappyStepModel,
    ToolModel,
)
from snappy_pipeline.wrapper_code import code_of, package_root
from snappy_pipeline.utils import listify
from snappy_pipeline.library_table import (
    TARGETED_LIBRARY_TYPES,
    _check_groups,
    _query_libraries,
    build_donor_dataframe,
    build_library_dataframe,
    colliding_libraries_message,
    entity_libraries,
    resolve_relationships,
    select_libraries,
)
from snappy_pipeline.reads import ReadGroup
from snappy_pipeline.workflows.abstract.protocol import (
    DATA_SETS,
    DataSignature,
    field_contract,
    select_signature,
)

if typing.TYPE_CHECKING:
    from snappy_pipeline.orchestration import Project
from snappy_pipeline.resource_usage import ResourceUsage, mem_mb, runtime_minutes


Inputs: typing.TypeAlias = InputFiles | dict[str, Any]
Outputs: typing.TypeAlias = OutputFiles | dict[str, Any]
Threads: typing.TypeAlias = int | None
Attempt: typing.TypeAlias = int | None


# ---------------------------------------------------------------------------
# Upstream group_by resolution
# ---------------------------------------------------------------------------
#: ``library_requirement`` of the steps and tools that compare each tumor with its matched
#: normal (the ``matched_normal_lib`` relationship); ``BaseStep.library_warnings`` names the
#: libraries they leave out
MATCHED_NORMAL_REQUIREMENT = "kind == 'cancer' and matched_normal_lib != ''"


#: The ``default_relationships`` of the steps that compare each tumor with its matched normal: the
#: donor's normal DNA library as ``matched_normal_lib``
MATCHED_NORMAL_RELATIONSHIP = {
    "matched_normal_lib": RelationshipDefinition(
        via="donor_name",
        target="role == 'normal' and extraction_type == 'dna'",
    )
}


@dataclass(frozen=True)
class RuleOptions:
    """The directives of a generated wrapper rule beyond the standard ones (``BaseStep.wrapper_rules``)"""

    #: ``localrule: True``: the job runs on the submit node
    local: bool = False
    #: ``retries``: how often a failed job is restarted (network access)
    retries: int = 0
    #: ``cache: "omit-software"``: the output is cached between workflows; such a rule has no benchmark
    cache: bool = False
    #: ``wildcard_constraints`` of the rule
    constraints: dict[str, str] = field(default_factory=dict)
    #: The step Snakefile writes the rule itself (a container, a checkpoint, an output that is no dict)
    explicit: bool = False


@dataclass(frozen=True)
class ReportOutput:
    """An output that the Snakemake report shows (``snappy report``)."""

    #: Subcategory below the task's category
    subcategory: str | None = None
    #: Caption, an RST file in the ``report/`` directory next to the step part's module
    caption: str | None = None
    #: For a directory output: the HTML file to open (Snakemake's ``htmlindex``)
    htmlindex: str | None = None
    #: For a directory output: the files to show, as patterns relative to it (Snakemake's
    #: ``patterns``); their wildcards label the files
    patterns: tuple[str, ...] = ()


class BaseStepPart:
    """Base class for a part of a pipeline step"""

    name: str

    #: The actions available in the class.
    actions: tuple[str, ...]

    #: The wrapper of each action that runs one, as a path below ``snappy_wrappers/wrappers``;
    #: a property where the config picks it. The shared include ``common/wrapper_rules.rules``
    #: writes the rule of each of them (``BaseStep.wrapper_rules``) with ``wf.get_wrapper(...)``,
    #: and the hash of the wrapper's code is the ``code`` param of the rule (its rerun trigger).
    wrappers: dict[str, str] = {}

    #: The directives of the generated rules beyond the standard ones, by action (``RuleOptions``)
    rule_options: dict[str, RuleOptions] = {}

    #: The action whose command gets the ``extra_arguments`` of the tool section; ``None``: the
    #: only action
    extra_arguments_action: str | None = None

    #: Wildcard of the output entity (library, tumor library, ...) in the output paths
    entity_wildcard: str = "library_name"

    #: Output keys whose files are result files, for all actions or by action; ``get_result_files``
    #: derives the result files from them
    result_keys: tuple[str, ...] | dict[str, tuple[str, ...]] = ()

    #: Outputs that the Snakemake report shows, by action and output key
    report_outputs: dict[str, dict[str, ReportOutput]] = {}

    #: Output keys by action of scatter fragments (per-region or per-chunk pieces that a gather
    #: rule joins): always ``temp()``
    temp_outputs: dict[str, tuple[str, ...]] = {}

    #: Output keys by action of per-sample intermediates that only later rules of the same task
    #: read: ``temp()`` when the task's ``delete_intermediates`` names ``intermediates``. Result
    #: files, contract outputs and project-wide resources are neither.
    intermediate_outputs: dict[str, tuple[str, ...]] = {}

    #: Output keys by action of the task's BAM file and its index, which another task may process
    #: further: ``temp()`` and no result files when ``delete_intermediates`` names ``bam``
    bam_outputs: dict[str, tuple[str, ...]] = {}

    #: Default resource usage for actions that are not given in ``resource_usage``.
    default_resource_usage: ResourceUsage = ResourceUsage(
        threads=1,
        runtime="1h",
        mem="2GB",  # 1h
    )

    #: Configure resource usage here that should not use the default resource usage from
    #: ``default_resource_usage``.
    resource_usage: dict[str, ResourceUsage] = {}

    def __init__[P: BaseStep](self, parent: P):
        self.name = self.__class__.name
        self.parent: P = parent
        self.config = parent.config
        self.w_config = parent.w_config

    def _validate_action(self, action: str):
        """Validate provided action

        Checks that the provided ``action`` is listed in the valid class actions list.

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)

    def get_resource_usage(self, action: str) -> ResourceUsage:
        """Return the resource usage for the given action: ``resource_usage[action]`` or
        ``default_resource_usage``. A step part whose usage depends on its config defines either
        as a property."""
        self._validate_action(action)
        return self.resource_usage.get(action, self.default_resource_usage)

    @staticmethod
    def get_default_partition() -> str | None:
        """Helper that returns the default partition."""
        return os.getenv("SNAPPY_PIPELINE_PARTITION")

    def get_resource(
        self, action: str, resource_name: str
    ) -> Callable[[Wildcards, InputFiles, Threads, Attempt], Any]:
        """Return the amount of resources to be allocated for the given action.

        :param action: The action to return the resource requirement for.
        :param resource_name: The name to return the resource for.
        """
        if resource_name not in ("threads", "runtime", "mem", "partition", "tmpdir"):
            raise ValueError(f"Invalid resource name: {resource_name}")

        def _get_resource(
            wildcards: Wildcards = None,
            input: InputFiles = None,
            threads: Threads = None,
            attempt: Attempt = 1,
        ) -> Any:
            resource_usage = self.get_resource_usage(action)
            if resource_name == "tmpdir" and not resource_usage.tmpdir:
                return self.parent.get_tmpdir()
            if resource_name == "partition" and not resource_usage.partition:
                return self.get_default_partition()
            # '_get_resource' is the primary input function,
            # if resource is callable we need to call it
            resource = getattr(resource_usage, resource_name)
            if callable(resource):
                return resource(wildcards=wildcards, input=input, threads=threads, attempt=attempt)
            if resource_name in ("mem", "runtime"):
                return self._scaled(resource_name, resource_usage, input, attempt)
            return resource

        return _get_resource

    def _scaled(self, name: str, usage: ResourceUsage, input: InputFiles, attempt: int) -> str:
        """Return memory or runtime: base, plus a term per GB of input, grown with retries.

        The value grows by ``resources.retry_factor`` of the project config with each retry.
        Without growth or input term, the declared value is returned as it is. Partition limits
        are left to SLURM and the executor plugin.
        """
        base, per_gb = getattr(usage, name), getattr(usage, f"{name}_per_gb_input")
        factor = self.w_config.resources.retry_factor ** ((attempt or 1) - 1)
        if factor == 1 and not per_gb:
            return base
        parse = mem_mb if name == "mem" else runtime_minutes
        value = parse(base)
        if per_gb:
            value += parse(per_gb) * input.size_mb / 1024
        value *= factor
        return f"{math.ceil(value)}MB" if name == "mem" else f"{math.ceil(value)}m"

    def tool_section(self) -> ToolModel | None:
        """Return the config section of the step part's tool, ``config.<name>``"""
        section = getattr(self.config, self.name, None)
        return section if isinstance(section, ToolModel) else None

    def with_extra_arguments(
        self, action: str, params: Callable[..., dict[str, Any]]
    ) -> Callable[..., dict[str, Any]]:
        """Add the ``extra_arguments`` of the tool section to the params of the action that runs
        the tool's command, unless the params already have them"""
        section = self.tool_section()
        target = self.extra_arguments_action or (
            self.actions[0] if len(self.actions) == 1 else None
        )
        if section is None or action != target:
            return params

        @wraps(params)
        def with_extra(*args, **kwargs):
            return {"extra_arguments": section.extra_arguments} | params(*args, **kwargs)

        return with_extra

    def get_params(self, action: str) -> Callable[..., dict[str, Any]]:
        """Return the params function of the given action: ``self._get_params_<action>``, or
        ``_no_params`` for an action without one.

        The function takes ``wildcards`` (and optionally ``input``) and returns the dict that
        the wrapper reads as ``snakemake.params.args``.
        """
        self._validate_action(action)
        return getattr(self, f"_get_params_{action}", self._no_params)

    @staticmethod
    def _no_params(wildcards) -> dict[str, Any]:
        """The params of an action that has none of its own (``extra_arguments`` and ``code`` are added anyway)"""
        return {}

    def get_input_files(self, action: str) -> Callable[[Wildcards], Inputs]:
        """Return the input function of the given action: ``self._get_input_files_<action>``, or
        ``_no_inputs`` for an action without inputs.

        The function takes ``wildcards`` and returns a dict (named inputs, used with
        ``unpack()``), a list or a single path.
        """
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}", self._no_inputs)

    @staticmethod
    def _no_inputs(wildcards) -> list[str]:
        """The inputs of an action that has none"""
        return []

    def get_output_files(self, action: str) -> Outputs:
        """Return the output files of an action: ``self._get_output_files_<action>()``."""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def get_log_file(self, action: str) -> Outputs:
        """Return the log files of an action: ``self._get_log_file_<action>()``, or
        ``self._get_log_file(action)`` for a step part whose actions share one pattern."""
        self._validate_action(action)
        if method := getattr(self, f"_get_log_file_{action}", None):
            return method()
        return self._get_log_file(action)

    def get_result_files(self) -> list[str]:
        """Return the result files, ``work/`` paths that the task's results job links into
        ``output/``: the outputs of ``result_keys`` for each output entity; a static path (a file
        of the whole task, such as a panel of normals) once."""
        result = []
        for action in self.actions:
            keys = self.result_keys
            if isinstance(keys, dict):
                keys = keys.get(action, ())
            if Intermediate.bam in self.config.delete_intermediates:
                keys = set(keys) - set(self.bam_outputs.get(action, ()))
            if not keys:
                continue
            for key, paths in self.get_output_files(action).items():
                if key not in keys:
                    continue
                for path in paths if isinstance(paths, list) else [paths]:
                    template = re.sub(r"\{(\w+),[^}]*\}", r"{\1}", str(path))
                    if "{" not in template:
                        result.append(template)
                        continue
                    for entity in self.parent.output_entities:
                        result += expand(template, **self._result_wildcards(entity))
        return result

    def _result_wildcards(self, entity: str) -> dict[str, list[str]]:
        """Return the wildcard values of the result files of one output entity."""
        return {self.entity_wildcard: [entity]}


@dataclass(frozen=True)
class ResolvedDependency:
    """Normalized dependency resolution result for one ``depends_on`` field."""

    field_name: str
    task_name: str
    step_name: str
    workflow_cls: type["BaseStep"]
    signature: DataSignature | None
    expected_schema: type[pydantic.BaseModel] | None


@dataclass(frozen=True)
class WrapperRule:
    """One generated wrapper rule: its name, the step part and action it wires, and its options"""

    name: str
    part: str
    action: str
    options: RuleOptions


class WritePedigreeStepPart(BaseStepPart):
    """Write the PED file of a pedigree, whose family column is the pedigree's name (its ``related_group``).

    The rows come from the sample sheets through the params, so that a changed sex, parent or
    affected status rewrites the PED and reruns its readers.
    """

    #: Step name
    name = "write_pedigree"

    #: Class available actions
    actions = ("run",)

    #: The wildcard that names the pedigree in the output path
    entity_wildcard = "related_group"

    #: The column that names the individuals and their parents: their DNA libraries, the
    #: sample names of the VCFs
    individual_column = "library_name"

    def get_result_files(self):
        return []

    @listify
    def _get_input_files_run(self, wildcards):
        """Return the BAM files of the pedigree members."""
        for library_name in self._members(getattr(wildcards, self.entity_wildcard))["library_name"]:
            yield self.parent.get_upstream_paths("alignments", library_name=library_name).bam

    def _get_output_files_run(self):
        return f"work/write_pedigree.{{{self.entity_wildcard}}}/out/{{{self.entity_wildcard}}}.ped"

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        """The PED rows as biomedsheets wrote them: the family, the individual, its father and mother, sex (1 male, 2 female, 0 unknown) and affected status (2 affected, 1 unaffected, 0 unknown)"""
        entity = getattr(wildcards, self.entity_wildcard)
        sexes = {"male": "1", "female": "2"}
        rows = [
            [
                self._family(entity),
                self._individual(row["library_name"]),
                self._individual(row["father_name"]),
                self._individual(row["mother_name"]),
                sexes.get(str(row["sex"]).lower(), "0"),
                str(row["disease_state"]),
            ]
            for _, row in self._members(entity).iterrows()
        ]
        return {"pedigree": rows}

    def _members(self, entity: str):
        """Return the primary DNA library rows of the pedigree ``entity``; none in projects with cancer sample sheets only."""
        members = self.parent.pedigree_members(entity)
        df = self.parent.build_library_dataframe()
        return members.iloc[0:0] if df.empty or df["kind"].eq("cancer").all() else members

    def _family(self, entity: str) -> str:
        """The family column of the PED: the pedigree's name, ``FAM_`` and the index donor on a germline sheet without a family id"""
        return entity

    def _individual(self, library_name: str) -> str:
        """The PED name of the individual of a DNA library, ``0`` for none"""
        if library_name == "0" or self.individual_column == "library_name":
            return library_name
        return self.parent.library_row(library_name)[self.individual_column]

    @staticmethod
    def write(args: dict[str, Any], output: OutputFiles):
        """Write the PED rows of the params (the rule's ``run:`` block)"""
        with open(str(output), "w") as f:
            for row in args["pedigree"]:
                print("\t".join(row), file=f)


class WritePedigreeSampleNameStepPart(WritePedigreeStepPart):
    """Write the pedigree with the donor names (``P001`` instead of ``P001-N1-DNA1-WGS1``)

    Used by export-external workflows that operate on externally-provided data, so they
    have no upstream ``ngs_mapping`` dependency and no BAM files to order by.
    """

    #: Step name
    name = "write_pedigree_with_sample_name"

    individual_column = "donor_name"

    def _get_input_files_run(self, wildcards):
        """Return no input files: pedigree writing only needs sample-sheet data."""
        return []


class LinkOutStepPart(BaseStepPart):
    """Generically link out

    This is for output files that are created unconditionally, i.e., for output files where the
    output name is the same as for the work file.
    """

    name = "link_out"
    actions = ("run",)

    def __init__(self, parent, disable_patterns=None):
        super().__init__(parent)
        self.base_pattern_out = "output/{path}/{file}.{ext}"
        self.base_path_out = self.base_pattern_out.replace(",[^/]", "")
        self.base_path_in = "work/{path}/{file}.{ext}"
        #: Patterns for disabling linking out to.  This is useful/required when there is a
        #: specialized link out step part, e.g., for the case of alignment where realignment is
        #: performed or not, depending on the configuration.
        self.disable_patterns = list(disable_patterns or [])

    def get_output_files(self, action):
        """Return output file pattern"""
        assert action == "run", "Unsupported action"
        return self.base_pattern_out

    def _get_input_files_run(self, wildcards):
        """Return the work file to link to"""
        path = self.base_path_in.format(**wildcards)
        for pattern in self.disable_patterns:
            if fnmatch(path, pattern):
                raise ValueError("Blocking match...")
        return path

    def run_locally(self, action, wildcards):
        assert action == "run", "Unsupported action"
        task_prefix = self.parent.task_path_prefix()
        # Prepend the task prefix to the paths
        path_out = task_prefix + f"output/{wildcards.path}/{wildcards.file}.{wildcards.ext}"
        path_in = task_prefix + f"work/{wildcards.path}/{wildcards.file}.{wildcards.ext}"
        if not os.path.islink(path_out):
            target = os.path.relpath(path_in, start=os.path.dirname(path_out))
            os.symlink(target, path_out)


class BaseStep:
    """Base class for the pipeline steps

    Each pipeline step is a Snakemake workflow
    """

    #: Override with step name
    name: str

    #: DataSignatures a task of this step produces; the default of :meth:`task_produces`
    produces: list[DataSignature] = []

    #: Override with the Pydantic model class for configuration validation
    config_model_class: type[SnappyStepModel]

    #: The step part of each tool, by the ``tool`` value of the task config: the constructor
    #: registers the selected one after the ``helper_parts``; a step that registers its parts
    #: itself leaves it empty
    tool_parts: dict[str, type[BaseStepPart]] = {}

    #: Structural step parts registered for every tool, such as ``WritePedigreeStepPart``
    helper_parts: tuple[type[BaseStepPart], ...] = ()

    #: The libraries that the step can process at all, as a pandas query, such as RNA libraries
    #: for expression quantification; it applies in addition to ``library_selection``.
    library_requirement: str | None = None

    #: The libraries that a task takes without ``library_selection``, as a pandas query.
    #: ``None`` takes the default of each library's sheet kind (cancer: tumor DNA, germline:
    #: DNA), ``""`` takes all libraries.
    library_default: str | None = None

    #: Default relationships merged into ``build_library_dataframe()``.
    #: Config-level ``relationships`` override these.
    default_relationships: dict[str, RelationshipDefinition] = {}

    #: Whether the step produces results for the libraries of background data sets, too
    background_results: bool = False

    #: Whether the step uses the libraries of background data sets otherwise, such as cohort
    #: members of joint calling or the samples of a model; only such steps and those with
    #: ``background_results`` accept ``background_selection``
    uses_background: bool = False

    #: Whether ``get_output_paths`` needs the project's read discovery (``discovery`` keyword)
    needs_discovery: bool = False

    def __init_subclass__(cls, **kwargs):
        """Enforce that every ``BaseStep`` subclass declares a ``produces`` contract."""
        super().__init_subclass__(**kwargs)
        # Only check concrete leaf classes that inherit produces=[] from BaseStep
        own_produces = cls.__dict__.get("produces")
        if own_produces is None and not cls.produces:
            raise TypeError(
                f"'{cls.__name__}' must declare a 'produces' list for contract compliance."
            )

    @classmethod
    def task_library_requirement(cls, config: SnappyStepModel) -> str | None:
        """Return the ``library_requirement`` of a task of this step; override it when the
        requirement depends on the task's config, such as its tool."""
        return cls.library_requirement

    @classmethod
    def task_library_default(cls, config: SnappyStepModel) -> str | None:
        """Return the ``library_default`` of a task of this step; override it when the default
        depends on the task's config, such as its tool."""
        return cls.library_default

    @classmethod
    def default_config_yaml(cls) -> str:
        """Return the annotated default configuration of a task of this step: a ``tasks:``
        list with one entry, generated from ``config_model_class``."""
        config = cls.config_model_class.default_config_yaml_string()
        lines = "".join(
            " " * 6 + line if line.strip() else line for line in config.splitlines(keepends=True)
        )
        return f"tasks:\n  - step: {cls.name}\n    name: {cls.name}\n    config:\n{lines}"

    @classmethod
    def task_produces(
        cls, config: SnappyStepModel, upstream: Mapping[str, DataSignature]
    ) -> tuple[DataSignature, ...]:
        """Return the DataSignatures that one task of this step produces.

        ``upstream`` maps each set ``depends_on`` field to the upstream signature this task reads
        (see :func:`~snappy_pipeline.workflows.abstract.protocol.select_signature`).
        The default returns :attr:`produces`; override it when the output depends on the task's
        config or its inputs. ``load_project()`` checks every ``depends_on`` requirement against
        the result.
        """
        return tuple(cls.produces)

    @classmethod
    def output_keys(
        cls, config: SnappyStepModel, signature: "DataSignature | None" = None
    ) -> set[str]:
        """Return the keys of ``get_output_paths`` for ``signature``.

        ``load_project()`` checks them against the consumers' ``ExpectedPathSchema``.
        """
        return set(cls.get_output_paths(config, signature))

    @classmethod
    def get_output_paths(
        cls, config: SnappyStepModel, signature: "DataSignature | None" = None, **kwargs
    ) -> "dict[str, str]":
        """Return **local** output paths for *signature*.

        Override in concrete workflow classes to expose this step's outputs to downstream
        consumers via :meth:`get_upstream_paths`.

        :param config: The config of the task whose outputs are requested.
        :param signature: The
            :class:`~snappy_pipeline.workflows.abstract.protocol.DataSignature` requested by the
            consumer. ``load_project()`` has already checked that the task produces it;
            implementations only use it to choose between several outputs.
        :param kwargs: Caller-supplied identifiers (e.g. ``library_name``, ``sample_name``)
            used to parametrise the returned path strings.  Implementations may use
            ``{library_name}``-style format strings when the identifiers are omitted.
        :returns: A ``dict[str, str]`` mapping logical output key (e.g. ``"bam"``, ``"vcf"``) to
            a **local** path string (relative to the task directory, e.g.
            ``"work/{library_name}/out/{library_name}.bam"``).
        :raises NotImplementedError: When the concrete subclass has not overridden this method.
        """
        raise NotImplementedError(
            f"'{cls.__name__}' must implement 'get_output_paths' to act as an upstream provider."
        )

    @classmethod
    def default_group_by(cls, config: SnappyStepModel) -> str | None:
        """Return the granularity of a task of this step that sets no ``group_by``, such as
        ``"related"`` for joint calling; ``None`` takes the one that upstream tasks declare."""
        return None

    @classmethod
    def filter_ids(cls, config: SnappyStepModel) -> frozenset[str]:
        """Return the ``FILTER`` IDs that a task of this step writes under names other than its
        task name; ``load_project()`` rejects other tasks named like one of them."""
        return frozenset()

    def __init__(self, workflow: Workflow, project: "Project", task_name: str):
        self.step_name = self.__class__.name
        self.workflow = workflow
        self.project = project
        self.task = project.task(task_name)
        if self.task.step != self.step_name:
            raise ValueError(
                f"Task {task_name!r} uses step {self.task.step!r}, not {self.step_name!r}"
            )
        self.task_name = task_name
        #: Setup logger for the step
        self.logger = logging.getLogger(task_name)
        self.logger.setLevel(logging.INFO if getattr(workflow, "verbose", False) else logging.WARN)

        self.w_config = project.model
        self.config = project.task_configs[task_name]
        self.depends_on = getattr(self.config, "depends_on", None)

        self.sub_steps: dict[str, BaseStepPart] = {}
        self._library_dataframes = {}
        self._primary_dna_libraries = {}
        #: The data sets with their sample sheets, loaded once per run and shared by all tasks
        self.data_set_infos = project.data_set_infos

        if self.tool_parts:
            self.register_sub_step_classes(
                (*self.helper_parts, self.tool_parts[str(self.config.tool)])
            )
        self.check_config()

    @staticmethod
    def task_root(task_name: str) -> str:
        """Return canonical namespace root for one task."""
        return f"tasks/{task_name}"

    @classmethod
    def namespaced_path(cls, task_name: str, local_path: str) -> str:
        """Return ``local_path`` namespaced below ``tasks/<task_name>/``."""
        return os.path.join(cls.task_root(task_name), local_path).replace("\\", "/")

    def task_path_prefix(self) -> str:
        """Return canonical task prefix for this workflow instance."""
        if getattr(self, "task_name", ""):
            return f"{self.task_root(self.task_name)}/"
        return ""

    def get_task_config(self, name: str) -> SnappyStepModel:
        """Return a task's validated config.

        ``name`` is this task's step or task name (returns its own config) or a ``depends_on``
        field (returns the config of the upstream task set there).
        """
        if name in (self.step_name, self.task_name):
            return self.config
        upstream = self.project.dependencies[self.task_name].get(name)
        if upstream is None:
            fields = list(type(self.depends_on).model_fields) if self.depends_on else []
            raise ValueError(
                f"Task {self.task_name!r}: depends_on.{name} names no task; "
                f"depends_on fields: {fields}"
            )
        if isinstance(upstream, tuple):
            raise ValueError(
                f"Task {self.task_name!r}: depends_on.{name} names several tasks; "
                "read their configs through upstream_task_names()"
            )
        return self.project.task_configs[upstream]

    def upstream_task_names(self, field: str) -> list[str]:
        """Return the task names of a set ``depends_on`` field: one, or several for a list field."""
        upstream = self.project.dependencies[self.task_name].get(field)
        if upstream is None:
            raise ValueError(f"Task {self.task_name!r}: depends_on.{field} names no task")
        return list(upstream) if isinstance(upstream, tuple) else [upstream]

    def read_groups(self, library_name: str, source: str | None = None) -> list[ReadGroup]:
        """Return the FASTQ read groups of ``library_name``, ordered by the path of the left file.

        ``source`` defaults to ``depends_on.reads``:

        - ``data_sets``: the folder named in the sample sheet, below the search paths of the
          library's data set;
        - an ``external_data`` task: a folder named like the library, below its search paths;
        - a task that writes FASTQs (``adapter_trimming``): its ``work/<library>/out``, where
          the files keep the names and sub-directories of the task's own input.
        """
        source = self.depends_on.reads if source is None else source
        if self.reads_from_task(source):
            out_dir = self.namespaced_path(source, f"work/{library_name}/out")
            upstream = self.project.task_configs[source].depends_on.reads
            return [group.rebased(out_dir) for group in self.read_groups(library_name, upstream)]
        folders, folder_name = self._library_folders(library_name, source)
        discovery = self.project.read_discovery
        groups = discovery.find(
            folders.search_paths, folder_name, folders.search_patterns, folders.mixed_se_pe
        )
        if not groups:
            raise ConfigError(discovery.missing(folders.search_paths, folder_name))
        return groups

    def _library_folders(self, library_name: str, source: str) -> tuple[LibraryFolders, str]:
        """Return where the reads of ``library_name`` are, and the name of its folder: the
        library's data set and the sheet's ``folderName``, or an ``external_data`` task and the
        library name."""
        if source in ("", DATA_SETS):
            row = self.library_row(library_name)
            folder_name = row.get("folder_name")
            if not isinstance(folder_name, str) or not folder_name:
                raise ConfigError(
                    f"Library {library_name!r} of data set {row['data_set']!r} has no folderName in "
                    "the sample sheet"
                )
            return self.w_config.data_sets[row["data_set"]], folder_name
        return self.project.task_configs[source], library_name

    def reads_input(self, library_name: str) -> list[str]:
        """Return the Snakemake inputs that provide the reads of ``library_name``.

        These are the FASTQ files themselves, or the ``out/`` directory of the task that writes
        them.
        """
        source = self.depends_on.reads
        if self.reads_from_task(source):
            return [self.namespaced_path(source, f"work/{library_name}/out")]
        return [path for group in self.read_groups(library_name) for path in group.paths.values()]

    def reads_from_task(self, source: str) -> bool:
        """Whether ``source`` (a ``depends_on.reads`` value) is a task that writes FASTQs, such as
        ``adapter_trimming``, rather than the data sets or an ``external_data`` task"""
        return source not in ("", DATA_SETS) and self.project.task(source).step != "external_data"

    def build_library_dataframe(self, include_background: bool = False):
        """Return the unified library dataframe for this workflow, with its relationships.

        ``include_background`` adds the libraries of background data sets (``is_background``).
        The table without relationships is built once per run and shared by all tasks; each
        task adds its relationships once. Callers must not modify the result in place.
        """
        if include_background not in self._library_dataframes:
            df = self._with_relationships(
                type(self), self.config, self.library_table(include_background)
            )
            self._library_dataframes[include_background] = self._narrow_background(self.config, df)
        return self._library_dataframes[include_background]

    @staticmethod
    def _narrow_background(config: SnappyStepModel, df):
        """Return ``df`` with only the background libraries that ``background_selection`` takes."""
        selection = getattr(config, "background_selection", None)
        if not selection or not df["is_background"].any():
            return df
        kept = set(_query_libraries(df[df["is_background"]], selection).index)
        return df[~df["is_background"] | df.index.isin(kept)]

    def library_table(self, include_background: bool):
        """Return the library table without relationships, built once per run and shared by all
        tasks (``Project.library_tables``)."""
        tables = self.project.library_tables
        if include_background not in tables:
            tables[include_background] = build_library_dataframe(
                self.data_set_infos, include_background=include_background
            )
        return tables[include_background]

    @staticmethod
    def _with_relationships(cls: type["BaseStep"], config: SnappyStepModel, table):
        """Return ``table`` with the relationship columns of a task of step ``cls``."""
        default_relationships = getattr(cls, "default_relationships", {})
        config_relationships = getattr(config, "relationships", None) or {}
        relationships = {**default_relationships, **config_relationships}
        return resolve_relationships(table, relationships) if relationships else table

    def task_libraries(self, task_name: str):
        """Return the rows of the libraries that the task ``task_name`` takes, like
        :attr:`libraries` of its workflow object, but computed from its config, so that it works
        for upstream tasks without a workflow object, too (``--frozen``)."""
        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        if task_name == self.task_name:
            return self.libraries
        cls = WORKFLOW_REGISTRY[self.project.task(task_name).step]
        config = self.project.task_configs[task_name]
        df = self._with_relationships(cls, config, self.library_table(cls.background_results))
        df = self._narrow_background(config, df)
        return select_libraries(
            df,
            getattr(config, "library_selection", None),
            cls.task_library_default(config),
            cls.task_library_requirement(config),
        )

    @cached_property
    def effective_group_by(self) -> str:
        """Return ``"library"`` or ``"related"``: the task's ``group_by``, else its step's
        :meth:`default_group_by`, else the one its upstream tasks declare (so a filtration
        after joint calling writes per group of related libraries), else ``"library"``."""
        return self.project.group_by[self.task_name]

    @cached_property
    def libraries(self):
        """Return the rows of the libraries that the task takes: the library table, with the
        background data sets if the step produces results for them (``background_results``),
        restricted to the step's ``library_requirement`` and selected by ``library_selection`` or
        else by the step's ``library_default`` (:func:`select_libraries`)."""
        df = self.build_library_dataframe(include_background=self.background_results)
        selection = getattr(self.config, "library_selection", None)
        return select_libraries(
            df,
            selection,
            self.task_library_default(self.config),
            self.task_library_requirement(self.config),
        )

    @cached_property
    def _entity_libraries(self) -> dict[str, str]:
        """Entity -> its primary library, see :func:`entity_libraries`."""
        return entity_libraries(self.libraries, self.effective_group_by)

    @cached_property
    def output_entities(self) -> list[str]:
        """Output entity names for this workflow: one per group of related libraries with ``group_by == "related"``, otherwise each of the task's libraries."""
        return list(self._entity_libraries)

    def primary_library(self, entity: str) -> str:
        """Return the primary library of an entity of this task: the library itself with ``group_by: library``, else the ``is_primary`` member of the related group (the pedigree index, the primary tumor) or its first one.
        Raises ``ValueError`` for a name that is no entity, such as a library behind a task that writes one result per related group."""
        if entity not in self._entity_libraries:
            raise ValueError(
                f"{entity!r} is no entity of task {self.task_name!r}, which writes one result per "
                f"{self.effective_group_by}: {', '.join(self._entity_libraries)}"
            )
        return self._entity_libraries[entity]

    def member_libraries(self, entity: str) -> list[str]:
        """Return the libraries that an entity of this task stands for in a joint result: with ``group_by: related`` the pedigree members of its related group (one DNA library per donor, the index first), else the library itself.
        A group without pedigree members, such as a cancer donor without a normal, stands for its primary library (the primary tumor).
        Raises ``ValueError`` for a name that is no entity."""
        if self.effective_group_by != "related":
            return [entity]
        members = self.pedigree_members(entity)["library_name"].tolist()
        if not members:
            return [self.primary_library(entity)]
        index = self.pedigree_index(entity)
        return [index, *(name for name in members if name != index)] if index else members

    def upstream_entity(self, field: str, library_name: str) -> str:
        """Return the entity under which the task of ``depends_on.<field>`` writes the result of ``library_name``: its related group when that task writes one result per related group, else the library itself."""
        (upstream,) = self.upstream_task_names(field)
        if self.project.group_by[upstream] != "related":
            return library_name
        return self.library_row(library_name)["related_group"]

    def primary_dna_libraries(self, include_background: bool = False):
        """Return the library dataframe rows of the primary DNA library of each donor.

        This is the germline view of the sample sheets: a donor's first DNA library that is not
        from a tumor sample. The column ``pedigree`` is the name of the donor's pedigree, its
        ``related_group``. On a cancer sample sheet, each donor is a pedigree of its own. A donor's
        RNA libraries and further DNA libraries are in the group, but no pedigree members. The
        task's ``library_selection`` applies first, so it narrows the pedigrees and their members.
        """
        if include_background not in self._primary_dna_libraries:
            df = self.build_library_dataframe(include_background=include_background)
            if selection := getattr(self.config, "library_selection", None):
                df = _query_libraries(df, selection)
            _check_groups(df)
            dna = df[(df["extraction_type"] == "dna") & (df["role"] != "tumor")]
            dna = dna.drop_duplicates(["related_group", "donor_name"])
            dna = dna.assign(pedigree=dna["related_group"])
            self._primary_dna_libraries[include_background] = dna
        return self._primary_dna_libraries[include_background]

    def pedigrees(self, include_background: bool = False) -> dict[str, list[str]]:
        """Return the primary DNA libraries of each pedigree, keyed by the pedigree's name.

        Pedigrees and their members keep the order of the sample sheets.
        """
        dna = self.primary_dna_libraries(include_background)
        return {
            index: list(members["library_name"])
            for index, members in dna.groupby("pedigree", sort=False)
        }

    def pedigree_members(self, pedigree: str, include_background: bool = False):
        """Return the primary DNA library rows of the pedigree named ``pedigree``."""
        dna = self.primary_dna_libraries(include_background)
        return dna[dna["pedigree"] == pedigree]

    def pedigree_of(self, library_name: str, include_background: bool = False) -> str:
        """Return the pedigree (the key of :meth:`pedigrees`) of the primary DNA library ``library_name``.

        Raises ``KeyError`` for a library that is no primary DNA library.
        """
        dna = self.primary_dna_libraries(include_background)
        pedigree = dna.loc[dna["library_name"] == library_name, "pedigree"]
        if pedigree.empty:
            raise KeyError(library_name)
        return pedigree.iloc[0]

    def pedigree_index(self, pedigree: str, include_background: bool = False) -> str | None:
        """Return the index library of ``pedigree``: its ``is_primary`` primary DNA library (the index, a cancer donor's normal), ``None`` when the index has no DNA library."""
        members = self.pedigree_members(pedigree, include_background)
        index = members[members["is_primary"]]
        return None if index.empty else index.iloc[0]["library_name"]

    def pedigree_donors(self, pedigree: str, include_background: bool = False):
        """Return the donor table rows of the pedigree named ``pedigree``: every donor of the
        pedigree, sequenced or not, with its parents as donor names, its sex and its DNA
        library (``""`` without one); see :func:`build_donor_dataframe`."""
        tables = self.project.donor_tables
        if include_background not in tables:
            tables[include_background] = build_donor_dataframe(
                self.data_set_infos, include_background=include_background
            )
        donors = tables[include_background]
        return donors[donors["pedigree"] == pedigree]

    def children_with_both_parents(self) -> list[str]:
        """Return the primary DNA libraries whose father and mother both have one."""
        dna = self.primary_dna_libraries()
        return dna[(dna["father_name"] != "0") & (dna["mother_name"] != "0")][
            "library_name"
        ].tolist()

    @cached_property
    def library_kits(self) -> dict[str, str]:
        """Return library name -> enrichment kit (``library_kit``) of the targeted libraries
        (``WES``, ``Panel_seq``) that have one. Other libraries use ``libraryKit`` for other
        things, such as the long-read chemistry that minimap2 reads."""
        kits = {}
        for name in self._library_rows:
            row = self.library_row(name)
            kit = row.get("library_kit")
            if row.get("library_type") in TARGETED_LIBRARY_TYPES and isinstance(kit, str) and kit:
                kits[name] = kit
        return kits

    def reference_inputs(self) -> dict[str, str]:
        """The genome FASTA of ``depends_on.reference`` and its ``.fai`` and ``.dict``, as rule
        inputs, so that a missing sidecar fails when the DAG is built."""
        reference = self.get_upstream_paths("reference")
        return {
            "reference": reference.fasta,
            "reference_fai": reference.fai,
            "reference_dict": reference.sequence_dict,
        }

    def vcf_inputs(self, field: str) -> dict[str, str]:
        """The VCF of ``depends_on.<field>`` and its index, as rule inputs ``<field>`` and
        ``<field>_index``, so that a missing index fails when the DAG is built."""
        vcf = self.get_upstream_paths(field)
        return {field: vcf.vcf, f"{field}_index": vcf.index}

    def library_kit_paths(self, field: str, library_name: str):
        """Return the files of the library's kit from the ``external_data`` task with ``kits``
        in ``depends_on.<field>``, like :meth:`get_upstream_paths`."""
        return self.get_upstream_paths(field, library_kit=self.library_kits[library_name])

    def check_library_kits(
        self, field: str, library_names: typing.Iterable[str], require: bool = True
    ) -> None:
        """Raise ``ConfigError`` unless the kits task in ``depends_on.<field>`` has the kit of
        every library; the message names all such libraries at once. With ``require=False``,
        libraries without a kit pass."""
        task_name = self.resolve_dependency(field).task_name
        kits = self.project.task_configs[task_name].kits
        missing = [
            f"{name} ({self.library_kits.get(name, 'no kit')})"
            for name in library_names
            if self.library_kits.get(name) not in kits and (require or name in self.library_kits)
        ]
        if missing:
            raise ConfigError(
                f"depends_on.{field} task {task_name!r} has no kit for these libraries: "
                f"{', '.join(missing)}; its kits: {', '.join(sorted(kits))}"
            )

    def common_kit(self, field: str, library_names: typing.Iterable[str]) -> str:
        """Return the one kit of the libraries, checked with :meth:`check_library_kits`; raise
        ``ConfigError`` when they have several."""
        library_names = list(library_names)
        self.check_library_kits(field, library_names)
        kits = sorted({self.library_kits[name] for name in library_names})
        if len(kits) != 1:
            raise ConfigError(
                f"the libraries must have one kit, but they have the kits "
                f"{', '.join(kits) or '(none)'}; use one task per kit"
            )
        return kits[0]

    @cached_property
    def _library_rows(self) -> dict[str, list[dict[str, typing.Any]]]:
        """Return library name -> its rows in all data sets, background included."""
        rows: dict[str, list[dict[str, typing.Any]]] = {}
        for row in self.build_library_dataframe(include_background=True).to_dict("records"):
            rows.setdefault(row["library_name"], []).append(row)
        return rows

    @cached_property
    def _selected_libraries(self) -> set[tuple[str, str]]:
        """Return (data set, library name) of the libraries that ``library_selection`` keeps."""
        df = self.build_library_dataframe(include_background=True)
        if selection := getattr(self.config, "library_selection", None):
            df = df.query(selection)  # apply_library_selection would reject colliding names
        return set(zip(df["data_set"], df["library_name"]))

    def library_row(self, library_name: str) -> dict[str, typing.Any]:
        """Return the row of ``library_name``: its data set, folder, kit and the other sample
        sheet columns. If several data sets have a library of that name, the task's
        ``library_selection`` must keep exactly one of them; otherwise this raises
        ``ConfigError`` with the fix."""
        rows = self._library_rows.get(library_name)
        if not rows:
            raise ConfigError(f"library {library_name!r} is not among the task's libraries")
        if len(rows) > 1:
            kept = [
                row for row in rows if (row["data_set"], library_name) in self._selected_libraries
            ]
            if len(kept) != 1:
                colliding = {library_name: [row["data_set"] for row in kept or rows]}
                raise ConfigError(colliding_libraries_message(colliding))
            rows = kept
        return rows[0]

    def matched(self, library_name: str, relationship: str) -> str:
        """The library that a relationship (``matched_normal_lib``, ``matched_rna_lib``) gives for
        ``library_name``, or ``""`` for none"""
        return self.library_row(library_name).get(relationship) or ""

    def matched_normal(self, tumor_library: str) -> str:
        """The matched normal DNA library of a tumor library (``MATCHED_NORMAL_RELATIONSHIP``), or ``""``"""
        return self.matched(tumor_library, "matched_normal_lib")

    def tumor_normal_alignments(self, tumor_library: str) -> dict[str, str]:
        """The rule inputs of a tool that compares a tumor with its matched normal: ``tumor_bam``,
        ``tumor_bai``, ``normal_bam`` and ``normal_bai`` from ``depends_on.alignments``"""
        tumor = self.get_upstream_paths("alignments", library_name=tumor_library)
        normal = self.get_upstream_paths(
            "alignments", library_name=self.matched_normal(tumor_library)
        )
        return {
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
            "normal_bam": normal.bam,
            "normal_bai": normal.bai,
        }

    def library_warnings(self) -> list[str]:
        """Return warnings about the libraries of this task, which the main Snakefile prints once
        per run: the cancer libraries without a matched normal that a task with
        ``MATCHED_NORMAL_REQUIREMENT`` leaves out, and the libraries of an explicit
        ``library_selection`` that the step's ``library_requirement`` leaves out."""
        requirement = self.task_library_requirement(self.config)
        if not requirement:
            return []
        warnings = []
        df = self.build_library_dataframe(include_background=self.background_results)
        selection = getattr(self.config, "library_selection", None)
        taken = set(self.libraries["library_name"])
        if requirement == MATCHED_NORMAL_REQUIREMENT:
            default = self.task_library_default(self.config)
            cancer = select_libraries(df, selection, default, "kind == 'cancer'")
            if unpaired := sorted(set(cancer["library_name"]) - taken):
                warnings.append(
                    f"Task {self.task_name!r} leaves out the libraries without a matched normal: "
                    + ", ".join(unpaired)
                )
                taken |= set(unpaired)
        if selection:
            selected = set(select_libraries(df, selection)["library_name"])
            if dropped := sorted(selected - taken):
                warnings.append(
                    f"Task {self.task_name!r} cannot process these libraries of its "
                    f"library_selection ({requirement}): " + ", ".join(dropped)
                )
        return warnings

    def check_config(self):
        """Check the configuration against the sample sheets at setup; raise ``ConfigError`` on
        problems. Override in sub classes."""

    def register_sub_step_classes(self, classes: tuple[type[BaseStepPart], ...]):
        """Register the step part classes: one object per class in ``self.sub_steps``, by name"""
        for klass in classes:
            self.sub_steps[klass.name] = klass(self)

    def resolve_dependency(
        self, field_name: str, upstream_task: str | None = None
    ) -> ResolvedDependency:
        """Resolve one typed ``depends_on`` field into a normalized dependency object.

        For a field that takes a list of tasks, ``upstream_task`` names one of them."""
        if self.depends_on is None:
            raise AttributeError(
                f"Task '{self.task_name}' has no 'depends_on' configuration; "
                f"cannot resolve upstream field '{field_name}'."
            )

        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        if field_name not in type(self.depends_on).model_fields:
            raise ValueError(
                f"Task '{self.task_name}' has no depends_on field named '{field_name}'."
            )
        # load_project checked the value and the upstream task; this only looks them up
        target_task_name = self.project.dependencies[self.task_name].get(field_name)
        if not target_task_name:
            raise ValueError(f"Task '{self.task_name}': depends_on.{field_name} is empty or unset.")
        if isinstance(target_task_name, tuple):
            if upstream_task not in target_task_name:
                raise ValueError(
                    f"Task '{self.task_name}': depends_on.{field_name} names several tasks; "
                    f"pass one of {list(target_task_name)} as upstream_task, not {upstream_task!r}"
                )
            target_task_name = upstream_task
        elif upstream_task not in (None, target_task_name):
            raise ValueError(
                f"Task '{self.task_name}': depends_on.{field_name} names {target_task_name!r}, "
                f"not {upstream_task!r}"
            )
        step_name = self.project.task(target_task_name).step
        contract = field_contract(
            type(self.depends_on).model_fields[field_name], self.config.tool_tag()
        )
        return ResolvedDependency(
            field_name=field_name,
            task_name=target_task_name,
            step_name=step_name,
            workflow_cls=WORKFLOW_REGISTRY[step_name],
            signature=contract.signature,
            expected_schema=contract.schema,
        )

    def upstream_signature(
        self, field_name: str, upstream_task: str | None = None
    ) -> DataSignature:
        """The signature that ``depends_on.<field_name>`` reads from its upstream task, as
        ``load_project`` selected it from what that task produces"""
        selected = self.project.selected[self.task_name][field_name]
        if isinstance(selected, tuple):
            names = self.project.dependencies[self.task_name][field_name]
            return selected[names.index(upstream_task)]
        return selected

    def _check_upstream_library(self, upstream: str, library_name: str) -> None:
        """Raise ``ConfigError`` if the upstream task writes no result for ``library_name``: a
        related group that it does not write (it writes per library, or its libraries leave the
        group out), or a library that it takes from another data set than this task (files are
        found by library name, so the upstream file would be of the other library)."""
        from snappy_pipeline.orchestration import _TASK_INSTANCES

        # Only between tasks of the running project; with --frozen, upstream tasks are not loaded
        instance = _TASK_INSTANCES.get(upstream)
        if _TASK_INSTANCES.get(self.task_name) is not self or instance is None:
            return
        if library_name not in self._library_rows:
            # A related group, not a library
            if library_name in instance.output_entities:
                return
            if instance.effective_group_by != "related":
                raise ConfigError(
                    f"task {upstream!r} writes one result per library, but {library_name!r} is a "
                    f"related group; set group_by: library on task {self.task_name!r}, or "
                    f"group_by: related on {upstream!r}"
                )
            raise ConfigError(
                f"task {upstream!r} writes no result for the related group {library_name!r}: its "
                f"libraries leave the group out; align library_selection of {self.task_name!r} "
                f"and {upstream!r}"
            )
        if library_name not in instance._library_rows:
            return
        own, theirs = (
            self.library_row(library_name)["data_set"],
            instance.library_row(library_name)["data_set"],
        )
        if own != theirs:
            raise ConfigError(
                f"library {library_name!r} comes from data set {own!r}, but the upstream task "
                f"{upstream!r} takes the library of that name from {theirs!r}. Select the "
                "libraries of the same data set in both tasks with library_selection."
            )

    def get_upstream_paths(
        self,
        req_field_name: str,
        signature: DataSignature | None = None,
        upstream_task: str | None = None,
        **kwargs,
    ) -> "pydantic.BaseModel | dict[str, str]":
        """Resolve global output paths from an upstream task using the Consumer-Driven Contract.

        This is the broker method for the semantic retrieval pattern.  The upstream
        task's workflow class is located via :data:`~snappy_pipeline.workflow_registry.WORKFLOW_REGISTRY`
        and its :py:meth:`get_output_paths` classmethod is called to produce **local** paths, which
        are then namespaced with the upstream task name.

        The ``depends_on`` field referenced by *req_field_name* may carry ``typing.Annotated``
        metadata with:

        * A :class:`~snappy_pipeline.workflows.abstract.protocol.DataSignature` – forwarded to
          :py:meth:`get_output_paths`.
        * An optional ``pydantic.BaseModel`` subclass – used to wrap the returned paths into a
          typed object for dot-notation access.

        Example::

            # in a step part's _get_input_files_run method
            alignments: ExpectedAlignments = self.parent.get_upstream_paths(
                "alignments", library_name=wildcards.library_name
            )
            bam = alignments.bam
            bai = alignments.bai

        Arguments:
            req_field_name: The field name on ``self.config.depends_on`` that holds the upstream
                task reference. Must match a key in the ``depends_on`` Pydantic model.
            signature: Narrows the field's requirement for this lookup, e.g. ``-filtered`` to
                read unfiltered calls. The provider's ``get_output_paths`` receives the produced
                signature that :func:`~snappy_pipeline.workflows.abstract.protocol.select_signature`
                picks for the requirement.
            upstream_task: For a field that takes a list of tasks, the one to read (see
                :py:meth:`upstream_task_names`).
            kwargs: Forwarded verbatim to the upstream workflow's
                :py:meth:`get_output_paths` classmethod (e.g. ``library_name``, ``sample_name``).

        Returns:
            A ``pydantic.BaseModel`` instance when an expected-schema class is present in the
            ``Annotated`` metadata, otherwise a plain ``dict[str, str]``.

        Raises:
            AttributeError: If ``self.config`` has no ``depends_on`` attribute.
            ValueError: If the upstream task or its workflow class cannot be resolved.
        """
        dependency = self.resolve_dependency(req_field_name, upstream_task)
        if "library_name" in kwargs:
            self._check_upstream_library(dependency.task_name, kwargs["library_name"])

        # The provider gets the signature this consumer reads, so it can choose the file.
        if signature is None:
            selected = self.upstream_signature(req_field_name, upstream_task)
        else:
            produced = self.project.signatures[dependency.task_name]
            selected = select_signature(produced, signature)
            if selected is None:
                raise ConfigError(
                    f"depends_on.{req_field_name} task {dependency.task_name!r} produces no "
                    f"{signature}"
                )
        if dependency.workflow_cls.needs_discovery:
            kwargs["discovery"] = self.project.read_discovery
        local_paths = dependency.workflow_cls.get_output_paths(
            self.project.task_configs[dependency.task_name], selected, **kwargs
        )

        # Task-local paths get the upstream task's prefix; absolute ones (existing files) do not.
        def global_path(path: str) -> str:
            return path if os.path.isabs(path) else self.namespaced_path(dependency.task_name, path)

        global_paths = {
            k: global_path(v) if isinstance(v, str) else [global_path(p) for p in v]
            for k, v in local_paths.items()
        }

        if dependency.expected_schema is not None:
            # The consumer's schema names what it needs; other outputs of the provider are ignored.
            fields = dependency.expected_schema.model_fields
            return dependency.expected_schema(
                **{key: path for key, path in global_paths.items() if key in fields}
            )
        return global_paths

    def get_input_files(self, sub_step: str, action: str) -> Callable[..., Inputs]:
        """Return the input function for action of substep, for the rule's ``input:`` section

        Delegates to the sub step object's get_input_files function
        """
        return self._require_function(sub_step, "get_input_files", action)

    def get_output_files(self, sub_step: str, action: str) -> Outputs:
        """Return list of strings with output files/patterns

        Delegates to the sub step object's get_output_files function, flags the outputs in its
        ``report_outputs`` for the Snakemake report and marks its ``temp_outputs`` (and, if the
        task deletes them, its ``intermediate_outputs`` and ``bam_outputs``) as ``temp()``.
        """
        part = self._get_sub_step(sub_step)
        outputs = part.get_output_files(action)
        if specs := part.report_outputs.get(action):
            outputs = {
                key: self._report(part, value, specs[key]) if key in specs else value
                for key, value in outputs.items()
            }
        keys = set(part.temp_outputs.get(action, ()))
        if Intermediate.intermediates in self.config.delete_intermediates:
            keys |= set(part.intermediate_outputs.get(action, ()))
        if Intermediate.bam in self.config.delete_intermediates:
            keys |= set(part.bam_outputs.get(action, ()))
        if keys:
            outputs = {key: temp(value) if key in keys else value for key, value in outputs.items()}
        return outputs

    def _report(self, part: BaseStepPart, path: str, spec: ReportOutput):
        """Flag ``path`` for the report, in the task's category, labelled by its wildcards."""
        caption = None
        if spec.caption:
            caption = os.path.join(
                os.path.dirname(inspect.getfile(type(part))), "report", spec.caption
            )
        wildcards = re.findall(r"\{(\w+)\}", " ".join((str(path), *spec.patterns)))
        labels = {name: "{%s}" % name for name in wildcards}
        return report(
            path,
            caption=caption,
            category=self.task_name,
            subcategory=spec.subcategory,
            labels=labels or None,
            htmlindex=spec.htmlindex,
            patterns=list(spec.patterns),
        )

    def get_benchmark_file(self, sub_step: str, action: str) -> str:
        """Return the benchmark file of a rule: the path of its log below ``benchmarks/``.

        Snakemake writes wall time, memory (``max_rss``), I/O and CPU time of the job into this
        TSV file; it calibrates the resource requests.
        """
        log = self.get_log_file(sub_step, action)
        if not log:  # a rule of an inactive tool, which has no outputs and never runs
            return f"benchmarks/inactive/{sub_step}.{action}.tsv"
        path = str(log["log"] if isinstance(log, dict) else log).removeprefix("work/")
        path = re.sub(r"(^|/)log/", r"\1", path, count=1)
        return f"benchmarks/{path.removesuffix('.log')}.tsv"

    def get_params(self, sub_step: str, action: str) -> Callable[..., Any]:
        """Return the params function for action of substep, for the rule's ``params:`` section

        Delegates to the sub step object's get_params function and adds the ``extra_arguments``
        of the tool section
        """
        params = self._require_function(sub_step, "get_params", action)
        part = self._get_sub_step(sub_step)
        with_arguments = part.with_extra_arguments(action, params)
        code = code_of(part.wrappers.get(action))
        if code is None:
            return with_arguments

        @wraps(with_arguments)  # Snakemake reads the signature to pick the arguments
        def with_code(*args, **kwargs):
            # The hash of the wrapper's code, so that a changed wrapper, tool or base class reruns
            # the job (Snakemake's own code trigger does not see wrapper scripts)
            return {**with_arguments(*args, **kwargs), "code": code}

        return with_code

    def _require_function(self, sub_step: str, method: str, action: str) -> Callable[..., Any]:
        """Call ``method`` of the sub step and check that it returned a function of wildcards"""
        part = self._get_sub_step(sub_step)
        value = getattr(part, method)(action)
        if not callable(value):
            raise TypeError(
                f"{type(part).__name__}.{method}({action!r}) must return a function of "
                f"wildcards, not {type(value).__name__}"
            )
        return value

    def get_resource(self, sub_step: str, action: str, resource_name: str) -> Any:
        """Return the function of one resource of a rule (``threads:`` takes ``"threads"``)"""
        return self._get_sub_step(sub_step).get_resource(action, resource_name)

    def get_resources(self, sub_step: str, action: str) -> dict[str, Callable[..., Any]]:
        """Return the functions of a rule's ``resources:`` section: runtime, mem, partition and
        tmpdir, for ``resources: **wf.get_resources(...)``"""
        return {
            name: self.get_resource(sub_step, action, name)
            for name in ("runtime", "mem", "partition", "tmpdir")
        }

    @cached_property
    def tmpdir(self) -> str:
        """The temporary directory of the jobs (the ``tmpdir`` resource): ``resources.tmpdir`` of the
        project config, else ``$TMPDIR``, else the system's, as for Snakemake. It is created once."""
        tmpdir = os.path.expandvars(self.w_config.resources.tmpdir) or os.getenv("TMPDIR")
        tmpdir = tmpdir or tempfile.gettempdir()
        os.makedirs(tmpdir, exist_ok=True)
        return tmpdir

    def get_tmpdir(self) -> str:
        """See ``tmpdir``; used by ``get_resource(part, action, "tmpdir")``"""
        return self.tmpdir

    def get_log_file(self, sub_step: str, action: str) -> Outputs:
        """Return the log files of a rule, for ``log: **wf.get_log_file(...)``"""
        return self._get_sub_step(sub_step).get_log_file(action)

    def wrapper_rules(self, cached: bool = False) -> list[WrapperRule]:
        """Return the rules that the shared include ``common/wrapper_rules.rules`` writes: one per
        action with a wrapper of every registered step part, named ``<step>_<part>_<action>``,
        unless the part marks it ``explicit``; ``cached`` selects the rules cached between
        workflows, which have no benchmark."""
        rules = []
        for part in self.sub_steps.values():
            for action in part.actions:
                if action not in part.wrappers:
                    continue
                options = part.rule_options.get(action, RuleOptions())
                if options.explicit or options.cache != cached:
                    continue
                name = f"{self.step_name}_{part.name}_{action}"
                rules.append(WrapperRule(name, part.name, action, options))
        return rules

    def get_wrapper(self, sub_step: str, action: str) -> str:
        """Return the wrapper of a rule as a ``file://`` URL: the step part's ``wrappers[action]``"""
        part = self._get_sub_step(sub_step)
        part._validate_action(action)
        if action not in part.wrappers:
            raise ValueError(f"{type(part).__name__}.wrappers names no wrapper for {action!r}")
        return self.wrapper_path(part.wrappers[action])

    def get_result_files(self) -> list[str]:
        """Return the result files (``work/`` paths) of all step parts."""
        return [path for part in self.sub_steps.values() for path in part.get_result_files()]

    def _get_sub_step(self, sub_step: str) -> BaseStepPart:
        if sub_step in self.sub_steps:
            return self.sub_steps[sub_step]
        else:
            raise ValueError(
                'Could not find sub step "{}" in workflow step "{}"'.format(sub_step, self.name)
            )  # pragma: no cover

    @staticmethod
    def wrapper_path(path: str) -> str:
        """The ``file://`` URL of a wrapper directory below ``snappy_wrappers/wrappers``"""
        return (package_root() / "wrappers" / path).as_uri()
