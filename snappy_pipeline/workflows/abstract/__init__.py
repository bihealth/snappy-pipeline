# -*- coding: utf-8 -*-
"""Base classes for the actual pipeline steps"""

import contextlib
import datetime
import logging
import os
import os.path
import inspect
import math
import re
import sys
import tempfile
import typing
from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass
from fnmatch import fnmatch
from functools import cached_property, lru_cache, wraps
from typing import Any, Callable

import pydantic
from biomedsheets import io_tsv
from biomedsheets.io import SheetBuilder, json_loads_ordered
from biomedsheets.io_tsv.base import LIBRARY_TYPES, SheetIOException
from biomedsheets.naming import NAMING_SCHEMES, name_generator_for_scheme
from biomedsheets.ref_resolver import RefResolver
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
from snappy_pipeline.reads import ReadGroup
from snappy_pipeline.workflows.abstract.protocol import (
    DATA_SETS,
    DataSignature,
    ExpectedPathSchema,
    select_signature,
)

if typing.TYPE_CHECKING:
    from snappy_pipeline.orchestration import Project
from snappy_wrappers.resource_usage import ResourceUsage, mem_mb, runtime_minutes


@contextlib.contextmanager
def modified_environ(*remove, **update):
    """
    Temporarily updates the ``os.environ`` dictionary in-place.

    The ``os.environ`` dictionary is updated in-place so that the modification
    is sure to work in all situations.

    :param remove: Environment variables to remove.
    :param update: Dictionary of environment variables and values to add/update.

       Source: https://stackoverflow.com/a/34333710/84349
    """
    env = os.environ
    update = update or {}
    remove = remove or []

    # List of environment variables being updated or removed.
    stomped = (set(update.keys()) | set(remove)) & set(env.keys())
    # Environment variables and values to restore on exit.
    update_after = {k: env[k] for k in stomped}
    # Environment variables and values to remove on exit.
    remove_after = frozenset(k for k in update if k not in env)

    try:
        env.update(update)
        [env.pop(k, None) for k in remove]
        yield
    finally:
        env.update(update_after)
        [env.pop(k) for k in remove_after]


Inputs: typing.TypeAlias = InputFiles | dict[str, Any]
Outputs: typing.TypeAlias = OutputFiles | dict[str, Any]
Threads: typing.TypeAlias = int | None
Attempt: typing.TypeAlias = int | None


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

    #: Default resource usage for actions that are not given in ``resource_usage``.
    default_resource_usage: ResourceUsage = ResourceUsage(
        threads=1,
        runtime="1h",
        mem="2GB",  # 1h
    )

    #: Configure resource usage here that should not use the default resource usage from
    #: ``default_resource_usage``.
    resource_usage: dict[str, ResourceUsage] = {}

    #: Outputs that the Snakemake report shows, by action and output key
    report_outputs: dict[str, dict[str, ReportOutput]] = {}

    #: Output keys whose files are result files, for all actions or by action; ``get_result_files``
    #: derives the result files from them
    result_keys: tuple[str, ...] | dict[str, tuple[str, ...]] = ()

    #: Wildcard of the output entity (library, tumor library, ...) in the output paths
    entity_wildcard: str = "library_name"

    #: The wrapper of each action that runs one, as a path below ``snappy_wrappers/wrappers``;
    #: a property where the config picks it. The rule's ``wrapper:`` is ``wf.get_wrapper(...)``,
    #: and the hash of the wrapper's code is the ``code`` param of the rule (its rerun trigger).
    wrappers: dict[str, str] = {}

    #: The action whose command gets the ``extra_arguments`` of the tool section; ``None``: the
    #: only action
    extra_arguments_action: str | None = None

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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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
            resource_usage = self.get_resource_usage(
                action, wildcards=wildcards, input=input, threads=threads, attempt=attempt
            )
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
        """Return the input function of the given action: ``self._get_input_files_<action>``

        The function takes ``wildcards`` and returns a dict (named inputs, used with
        ``unpack()``), a list or a single path.
        """
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

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
        """Return the result files below ``output/``: the outputs of ``result_keys`` for
        each output entity, at their mirror path."""
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
                    unconstrained = re.sub(r"\{(\w+),[^}]*\}", r"{\1}", str(path))
                    template = "output/" + unconstrained.removeprefix("work/")
                    for entity in self.parent.output_entities:
                        result += expand(template, **self._result_wildcards(entity))
        return result

    def _result_wildcards(self, entity: str) -> dict[str, list[str]]:
        """Return the wildcard values of the result files of one output entity."""
        return {self.entity_wildcard: [entity]}


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

    def _members(self, entity: str):
        """Return the primary DNA library rows of the pedigree ``entity``; none in projects with cancer sample sheets only."""
        members = self.parent.pedigree_members(entity)
        df = self.parent.build_library_dataframe()
        return members.iloc[0:0] if df.empty or df["kind"].eq("cancer").all() else members

    def _family(self, entity: str) -> str:
        """The family column of the PED: the pedigree's name, ``FAM_`` and the index donor on a germline sheet without a family id"""
        return entity

    @listify
    def _get_input_files_run(self, wildcards):
        """Return the BAM files of the pedigree members."""
        for library_name in self._members(getattr(wildcards, self.entity_wildcard))["library_name"]:
            yield self.parent.get_upstream_paths("alignments", library_name=library_name).bam

    def _get_output_files_run(self):
        return f"work/write_pedigree.{{{self.entity_wildcard}}}/out/{{{self.entity_wildcard}}}.ped"

    def get_result_files(self):
        return []

    def _individual(self, library_name: str) -> str:
        """The PED name of the individual of a DNA library, ``0`` for none"""
        if library_name == "0" or self.individual_column == "library_name":
            return library_name
        return self.parent.library_row(library_name)[self.individual_column]

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

    def _get_input_files_run(self, wildcards):
        """Return the work file to link to"""
        path = self.base_path_in.format(**wildcards)
        for pattern in self.disable_patterns:
            if fnmatch(path, pattern):
                raise ValueError("Blocking match...")
        return path

    def get_output_files(self, action):
        """Return output file pattern"""
        assert action == "run", "Unsupported action"
        return self.base_pattern_out

    def run_locally(self, action, wildcards):
        assert action == "run", "Unsupported action"
        task_prefix = self.parent.task_path_prefix()
        # Prepend the task prefix to the paths
        path_out = task_prefix + f"output/{wildcards.path}/{wildcards.file}.{wildcards.ext}"
        path_in = task_prefix + f"work/{wildcards.path}/{wildcards.file}.{wildcards.ext}"
        if not os.path.islink(path_out):
            target = os.path.relpath(path_in, start=os.path.dirname(path_out))
            os.symlink(target, path_out)


@lru_cache()
def _cached_read_cancer_tsv_sheet(path_abs, path_rel, naming_scheme):
    """Cached reading of cancer TSV sheets."""
    with open(path_abs, "rt") as f:
        return io_tsv.read_cancer_tsv_sheet(f, path_rel, naming_scheme)


@lru_cache()
def _cached_read_germline_tsv_sheet(path_abs, path_rel, naming_scheme):
    """Cached reading of germline TSV sheets."""
    with open(path_abs, "rt") as f:
        return io_tsv.read_germline_tsv_sheet(f, path_rel, naming_scheme)


@lru_cache()
def _cached_read_generic_tsv_sheet(path_abs, path_rel, naming_scheme):
    """Cached reading of generic TSV sheets."""
    with open(path_abs, "rt") as f:
        return io_tsv.read_generic_tsv_sheet(f, path_rel, naming_scheme)


class DataSetInfo:
    """Information on a DataSet"""

    def __init__(
        self,
        name,
        sheet_path,
        sheet_type,
        is_background,
        naming_scheme,
        sodar_uuid,
        sodar_title,
        pedigree_field=None,
    ):
        """Constructor.

        :param name: Name of the data set.
        :type name: str

        :param sheet_path: Path to sheet file that should be loaded.
        :type sheet_path: str

        :param sheet_type: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will
        attempt to load from sheet.
        :type sheet_type: str

        :param is_background: Whether the data set info is to be used only for background.
        :type is_background: bool

        :param naming_scheme: Selected naming schema: either 'secondary_id_pk' or
        'only_secondary_id'.
        :type naming_scheme: str

        :param sodar_uuid: The UUID of the corresponding SODAR project.
        :type sodar_uuid: str

        :param sodar_title: The title of the project in SODAR [optional].
        :type sodar_title: str

        :param pedigree_field: Custom field from sample sheet used to define pedigree, e.g.,
        'familyId'. If none defined, it will set pedigree based on sample sheet 'row'.
        Default: None.
        :type pedigree_field: str
        """
        #: Name of the data set
        self.name = name
        #: Path to the sheet file, for loading
        self.sheet_path = sheet_path
        #: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will attempt to load
        # from sheet.
        self.sheet_type = sheet_type
        #: Whether the data set info is to be used only for background
        self.is_background = is_background
        #: Selected naming schema
        if naming_scheme not in NAMING_SCHEMES:
            raise ValueError("Invalid naming scheme: {}".format(naming_scheme))  # pragma: no cover
        self.naming_scheme = naming_scheme
        #: The BioMed SampleSheet
        self.sheet = self._load_sheet()
        #: The UUID of the corresponding SODAR project.
        self.sodar_uuid = sodar_uuid
        #: The (optional) title of the project in SODAR.
        self.sodar_title = sodar_title
        #: The (optional) custom field used to define pedigree
        self.pedigree_field_kwargs = None
        if pedigree_field:
            self.pedigree_field_kwargs = {"join_by_field": pedigree_field}

    def _load_sheet(self):
        # ``load_project()`` resolved the path against the project directory
        fname = self.sheet_path
        if not os.path.exists(fname):
            raise ConfigError(f"Could not find sample sheet file {fname}")
        if self.sheet_path.endswith(".tsv"):
            readers = {
                "matched_cancer": _cached_read_cancer_tsv_sheet,
                "germline_variants": _cached_read_germline_tsv_sheet,
                "generic": _cached_read_generic_tsv_sheet,
            }
            if self.sheet_type not in readers:
                raise ConfigError("Invalid sheet type {}".format(self.sheet_type))
            try:
                sheet = readers[self.sheet_type](fname, self.sheet_path, self.naming_scheme)
            except SheetIOException as e:
                raise ConfigError(f"Sample sheet {fname}: {e}") from None
            # biomedsheets skips rows without libraryType silently (germline and generic sheets)
            if not any(
                test_sample.ngs_libraries
                for donor in sheet.bio_entities.values()
                for sample in donor.bio_samples.values()
                for test_sample in sample.test_samples.values()
            ):
                raise ConfigError(
                    f"Sample sheet {fname}: no row has a libraryType ({', '.join(LIBRARY_TYPES)}), "
                    "so there are no libraries. A row without libraryType is a person without "
                    "sequencing data, such as a parent who was not sequenced."
                )
            return self._set_is_background(sheet, self.is_background)
        elif self.sheet_path.endswith(".json"):
            with open(fname, "rt") as f:
                sheet_json = json_loads_ordered(f.read())
            resolver = RefResolver(
                lookup_paths=os.path.dirname(os.path.abspath(self.sheet_path)),
                dict_class=OrderedDict,
            )
            resolved_json = resolver.resolve("file://" + self.sheet_path, sheet_json)
            return self._set_is_background(
                SheetBuilder(resolved_json).run(
                    name_generator=name_generator_for_scheme(self.naming_scheme)
                ),
                self.is_background,
            )
        else:
            raise ValueError(  # pragma: no cover
                "Invalid sheet file type of {}".format(self.sheet_path)
            )

    @classmethod
    def _set_is_background(cls, sheet, flag):
        """Override "is_background" flag"""
        # TODO: check whether is already there and fail if not compatible
        sheet.json_data["extraInfoDefs"]["is_background"] = {
            "type": "boolean",
            "default": False,
        }
        sheet.extra_infos["is_background"] = flag
        return sheet


@dataclass(frozen=True)
class ResolvedDependency:
    """Normalized dependency resolution result for one ``depends_on`` field."""

    field_name: str
    task_name: str
    step_name: str
    workflow_cls: type["BaseStep"]
    signature: DataSignature | None
    expected_schema: type[pydantic.BaseModel] | None


class BaseStep:
    """Base class for the pipeline steps

    Each pipeline step is a Snakemake workflow
    """

    #: Override with step name
    name: str

    #: DataSignatures a task of this step produces; the default of :meth:`task_produces`
    produces: list[DataSignature] = []

    #: Whether ``get_output_paths`` needs the project's read discovery (``discovery`` keyword)
    needs_discovery: bool = False

    #: Default relationships merged into ``build_library_dataframe()``.
    #: Config-level ``relationships`` override these.
    default_relationships: dict[str, RelationshipDefinition] = {}

    def __init_subclass__(cls, **kwargs):
        """Enforce that every ``BaseStep`` subclass declares a ``produces`` contract."""
        super().__init_subclass__(**kwargs)
        # Only check concrete leaf classes that inherit produces=[] from BaseStep
        own_produces = cls.__dict__.get("produces")
        if own_produces is None and not cls.produces:
            raise TypeError(
                f"'{cls.__name__}' must declare a 'produces' list for contract compliance."
            )

    #: The libraries that the step can process at all, as a pandas query, such as RNA libraries
    #: for expression quantification; it applies in addition to ``library_selection``.
    library_requirement: str | None = None

    #: The libraries that a task takes without ``library_selection``, as a pandas query.
    #: ``None`` takes the default of each library's sheet kind (cancer: tumor DNA, germline:
    #: DNA), ``""`` takes all libraries.
    library_default: str | None = None

    #: Whether the step produces results for the libraries of background data sets, too
    background_results: bool = False

    #: Whether the step uses the libraries of background data sets otherwise, such as cohort
    #: members of joint calling or the samples of a model; only such steps and those with
    #: ``background_results`` accept ``background_selection``
    uses_background: bool = False

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

    #: Override with the Pydantic model class for configuration validation
    config_model_class: type[SnappyStepModel]

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
        self.config_paths = project.config_paths
        self.project_dir = project.project_dir

        self.sub_steps: dict[str, BaseStepPart] = {}
        self._library_dataframes = {}
        self._primary_dna_libraries = {}
        self._donor_tables = {}
        self.data_set_infos = list(self._load_data_set_infos())
        self.sheets = [info.sheet for info in self.data_set_infos]

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
        if source not in ("", DATA_SETS) and self.project.task(source).step != "external_data":
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
        if source not in ("", DATA_SETS) and self.project.task(source).step != "external_data":
            return [self.namespaced_path(source, f"work/{library_name}/out")]
        return [path for group in self.read_groups(library_name) for path in group.paths.values()]

    def build_library_dataframe(self, include_background: bool = False):
        """Return the unified library dataframe for this workflow, with its relationships.

        ``include_background`` adds the libraries of background data sets (``is_background``).
        The table without relationships is built once per run and shared by all tasks; each
        task adds its relationships once. Callers must not modify the result in place.
        """
        if include_background not in self._library_dataframes:
            df = self._with_relationships(
                type(self), self.config, self._library_table(include_background)
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

    def _library_table(self, include_background: bool):
        """Return the library table without relationships, built once per run."""
        tables = self.project.library_tables
        if include_background not in tables:
            tables[include_background] = build_library_dataframe(
                self.data_set_infos, self.sheets, include_background=include_background
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
        df = self._with_relationships(cls, config, self._library_table(cls.background_results))
        df = self._narrow_background(config, df)
        return select_libraries(
            df,
            getattr(config, "library_selection", None),
            cls.task_library_default(config),
            cls.task_library_requirement(config),
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

    @cached_property
    def effective_group_by(self) -> str:
        """Return ``"library"`` or ``"related"``: the task's ``group_by``, else its step's
        :meth:`default_group_by`, else the one its upstream tasks declare (so a filtration
        after joint calling writes per group of related libraries), else ``"library"``."""
        own = getattr(self.config, "group_by", None)
        if own is not None:
            return own
        return _declared_group_by(self.project, self.task_name, frozenset()) or "library"

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
        if _declared_group_by(self.project, upstream, frozenset()) != "related":
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

    def pedigree_donors(self, pedigree: str, include_background: bool = True):
        """Return the donor table rows of the pedigree named ``pedigree``: every donor of the
        pedigree, sequenced or not, with its parents as donor names, its sex and its DNA
        library (``""`` without one); see :func:`build_donor_dataframe`."""
        if include_background not in self._donor_tables:
            self._donor_tables[include_background] = build_donor_dataframe(
                self.data_set_infos, self.sheets, include_background=include_background
            )
        donors = self._donor_tables[include_background]
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

    def library_warnings(self) -> list[str]:
        """Return warnings about the libraries of this task, which the main Snakefile prints once
        per run: by default, the cancer libraries without a matched normal that a task with
        ``MATCHED_NORMAL_REQUIREMENT`` leaves out."""
        if self.task_library_requirement(self.config) != MATCHED_NORMAL_REQUIREMENT:
            return []
        df = self.build_library_dataframe(include_background=self.background_results)
        selection = getattr(self.config, "library_selection", None)
        default = self.task_library_default(self.config)
        cancer = select_libraries(df, selection, default, "kind == 'cancer'")
        if cancer.empty:
            return []
        unpaired = sorted(set(cancer["library_name"]) - set(self.libraries["library_name"]))
        if not unpaired:
            return []
        return [
            f"Task {self.task_name!r} leaves out the libraries without a matched normal: "
            + ", ".join(unpaired)
        ]

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

        model_fields = type(self.depends_on).model_fields
        field_info = model_fields.get(field_name)
        if field_info is None:
            raise ValueError(
                f"Task '{self.task_name}' has no depends_on field named '{field_name}'."
            )

        target_task_name = getattr(self.depends_on, field_name)
        if not target_task_name:
            raise ValueError(f"Task '{self.task_name}': depends_on.{field_name} is empty or unset.")
        if isinstance(target_task_name, list):
            if upstream_task not in target_task_name:
                raise ValueError(
                    f"Task '{self.task_name}': depends_on.{field_name} names several tasks; "
                    f"pass one of {target_task_name} as upstream_task, not {upstream_task!r}"
                )
            target_task_name = upstream_task
        elif upstream_task not in (None, target_task_name):
            raise ValueError(
                f"Task '{self.task_name}': depends_on.{field_name} names {target_task_name!r}, "
                f"not {upstream_task!r}"
            )

        target_task = next((t for t in self.w_config.tasks if t.name == target_task_name), None)
        if target_task is None:
            raise ValueError(
                f"Task '{self.task_name}': upstream task '{target_task_name}' from "
                f"depends_on.{field_name} was not found in configuration."
            )
        step_name = target_task.step

        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        workflow_cls = WORKFLOW_REGISTRY.get(step_name)
        if workflow_cls is None:
            raise ValueError(
                f"Workflow class for step '{step_name}' (from depends_on.{field_name}="
                f"'{target_task_name}') not found in WORKFLOW_REGISTRY."
            )

        required_sig: DataSignature | None = None
        expected_schema: type[pydantic.BaseModel] | None = None
        tool = self.config.tool_tag()
        for meta in field_info.metadata:
            if isinstance(meta, DataSignature):
                required_sig = meta.for_tool(tool)
            elif isinstance(meta, ExpectedPathSchema):
                expected_schema = meta.for_tool(tool)
            elif isinstance(meta, type) and issubclass(meta, pydantic.BaseModel):
                expected_schema = meta

        return ResolvedDependency(
            field_name=field_name,
            task_name=target_task_name,
            step_name=step_name,
            workflow_cls=workflow_cls,
            signature=required_sig,
            expected_schema=expected_schema,
        )

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
        required = signature or dependency.signature
        produced = self.project.signatures[dependency.task_name]
        selected = select_signature(produced, required)
        if selected is None:
            raise ConfigError(
                f"depends_on.{req_field_name} task {dependency.task_name!r} produces no {required}"
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

    def get_tmpdir(self) -> str:
        """Return temporary directory.

        To be used directly or via get_resource("step", "action", "tmpdir")

        1. Try to evaluate global_config/tmpdir. Interpret $-variables from environment.
           Provides the current date as $TODAY.
        2. If this fails, try to use environment variable TMPDIR.
        3. If this fails, use tempfile.gettempdir(), same as Snakemake default.
        """
        tmpdir = getattr(self.w_config, "global_config", {}).get("tmpdir", None)
        if tmpdir:
            with modified_environ(TODAY=datetime.date.today().strftime("%Y%m%d")):
                tmpdir = os.path.expandvars(tmpdir)
        if not tmpdir:
            tmpdir = os.getenv("TMPDIR")
        if not tmpdir:
            tmpdir = tempfile.gettempdir()
        # Force existence of temporary directory
        os.makedirs(tmpdir, exist_ok=True)
        return tmpdir

    def get_log_file(self, sub_step: str, action: str) -> Outputs:
        """Return the log files of a rule, for ``log: **wf.get_log_file(...)``"""
        return self._get_sub_step(sub_step).get_log_file(action)

    def get_wrapper(self, sub_step: str, action: str) -> str:
        """Return the wrapper of a rule as a ``file://`` URL: the step part's ``wrappers[action]``"""
        part = self._get_sub_step(sub_step)
        part._validate_action(action)
        if action not in part.wrappers:
            raise ValueError(f"{type(part).__name__}.wrappers names no wrapper for {action!r}")
        return self.wrapper_path(part.wrappers[action])

    def get_result_files(self) -> list[str]:
        """Return the result files below ``output/``: those of all step parts."""
        return [path for part in self.sub_steps.values() for path in part.get_result_files()]

    def _get_sub_step(self, sub_step: str) -> BaseStepPart:
        if sub_step in self.sub_steps:
            return self.sub_steps[sub_step]
        else:
            raise ValueError(
                'Could not find sub step "{}" in workflow step "{}"'.format(sub_step, self.name)
            )  # pragma: no cover

    def _load_data_set_infos(self) -> typing.Generator[DataSetInfo, None, None]:
        """Load BioMed Sample Sheets as given by configuration and yield them"""
        for name, data_set in self.w_config.data_sets.items():
            yield DataSetInfo(
                name,
                data_set.file,
                data_set.type,
                data_set.is_background,
                data_set.naming_scheme,
                data_set.sodar_uuid,
                data_set.sodar_title,
                data_set.pedigree_field,
            )

    @staticmethod
    def wrapper_path(path: str) -> str:
        """The ``file://`` URL of a wrapper directory below ``snappy_wrappers/wrappers``"""
        return (package_root() / "wrappers" / path).as_uri()


# ---------------------------------------------------------------------------
# Upstream group_by resolution
# ---------------------------------------------------------------------------


def _declared_group_by(project, task_name: str, visited: frozenset[str]) -> str | None:
    """Return the ``group_by`` that a task declares: its own, else its step's default, else the
    one that its upstream tasks declare, depth first in ``depends_on`` order. A task without any
    declares ``None``, so the final default ``"library"`` does not spread downstream."""
    from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

    config = project.task_configs[task_name]
    if (own := getattr(config, "group_by", None)) is not None:
        return own
    if default := WORKFLOW_REGISTRY[project.task(task_name).step].default_group_by(config):
        return default
    from snappy_pipeline.orchestration import upstream_tasks

    for upstream in upstream_tasks(project.dependencies.get(task_name, {})):
        if upstream not in visited and upstream != task_name:
            if declared := _declared_group_by(project, upstream, visited | {task_name}):
                return declared
    return None


# ---------------------------------------------------------------------------
# Library DataFrame helpers (used by variant_annotation, variant_filtration, …)
# ---------------------------------------------------------------------------

#: Default pandas query expressions applied when library_selection is None.
_LIBRARY_SELECTION_DEFAULTS: dict[str, str] = {
    "cancer": "role == 'tumor' and extraction_type == 'dna'",
    "germline": "extraction_type == 'dna'",
}

#: Library types whose ``libraryKit`` is an enrichment kit
TARGETED_LIBRARY_TYPES = ("WES", "Panel_seq")

#: ``library_requirement`` of the steps and tools that compare each tumor with its matched
#: normal (the ``matched_normal_lib`` relationship); ``BaseStep.library_warnings`` names the
#: libraries they leave out
MATCHED_NORMAL_REQUIREMENT = "kind == 'cancer' and matched_normal_lib != ''"

#: Standard columns of the library dataframe built by :func:`build_library_dataframe`.
_LIBRARY_COLUMNS = [
    "library_name",
    "extraction_type",
    "kind",
    "role",
    "is_primary",
    "sex",
    "donor_name",
    "sample_name",
    "related_group",
    "father_name",
    "mother_name",
    "disease_state",
    "tissue_type",
    "is_background",
    "data_set",
    "secondary_id",
]

#: Sample sheet keys that already feed standard columns and are not added again.
_SHEET_KEYS_IN_LIBRARY_COLUMNS = frozenset(
    (
        "extractionType",
        "isTumor",
        "sex",
        "isAffected",
        "fatherName",
        "motherName",
        "fatherPk",
        "motherPk",
    )
)


def _snake_case(name: str) -> str:
    """Return ``name`` in snake_case, e.g. ``libraryKit`` -> ``library_kit``."""
    name = re.sub(r"\W+", "_", name)
    name = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    return re.sub(r"_+", "_", name).strip("_").lower()


def _library_extra_columns(raw_sheet) -> dict[str, dict[str, Any]]:
    """Return library name -> extra columns from the sample sheet's ``extra_infos``.

    Each library gets the ``extra_infos`` of its bio entity, bio sample, test sample and
    itself, with keys converted to snake_case. Keys that already feed standard columns
    (:data:`_SHEET_KEYS_IN_LIBRARY_COLUMNS`) are skipped.

    Raises ``ConfigError`` when a key maps to a standard column, or when two keys map to
    the same column for one library with different values.
    """
    result: dict[str, dict[str, Any]] = {}
    for bio_entity in raw_sheet.bio_entities.values():
        for bio_sample in bio_entity.bio_samples.values():
            for test_sample in bio_sample.test_samples.values():
                for library in test_sample.ngs_libraries.values():
                    columns: dict[str, Any] = {}
                    sources: dict[str, str] = {}
                    levels = (
                        ("bio entity", bio_entity),
                        ("bio sample", bio_sample),
                        ("test sample", test_sample),
                        ("NGS library", library),
                    )
                    for level, entity in levels:
                        for key, value in (entity.extra_infos or {}).items():
                            if key in _SHEET_KEYS_IN_LIBRARY_COLUMNS:
                                continue
                            column = _snake_case(key)
                            source = f"{key!r} of {level} {entity.name!r}"
                            if column in _LIBRARY_COLUMNS:
                                raise ConfigError(
                                    f"Sample sheet column {source} maps to {column!r}, "
                                    "which is a standard library column."
                                )
                            if column in columns and columns[column] != value:
                                raise ConfigError(
                                    f"Sample sheet columns {sources[column]} and {source} both "
                                    f"map to {column!r} for library {library.name!r}, with "
                                    "different values."
                                )
                            columns.setdefault(column, value)
                            sources.setdefault(column, source)
                    result[library.name] = columns
    return result


def build_library_dataframe(
    data_set_infos,
    sheets,
    relationships=None,
    include_background=False,
):
    """Build a tidy pandas DataFrame with one row per NGS library.

    The biomedsheets views (CancerCaseSheet for cancer sheets; GenericSampleSheet for the
    libraries and GermlineCaseSheet for the pedigrees of the other sheets) are used directly,
    so that information already computed by biomedsheets (primary pairs, pedigree roles,
    affected status, sex, …) is re-used rather than re-derived.

    Parameters
    ----------
    data_set_infos:
        Iterable of :class:`DataSetInfo` objects (``self.data_set_infos``).
    sheets:
        Corresponding raw biomedsheets ``sheet`` objects (``self.sheets``).

    Returns
    -------
    pandas.DataFrame
        One row per NGS library.  Columns:

        * ``library_name`` — NGS library identifier
        * ``extraction_type`` — ``"dna"``, ``"rna"``, … (always lower-cased)
        * ``kind`` — ``"cancer"`` or ``"germline"``
        * ``role`` — ``"tumor"`` / ``"normal"`` / ``"index"`` / ``"father"``
          / ``"mother"`` / ``"affected"`` / ``"unaffected"``
        * ``is_primary`` — ``True`` for the primary tumor (cancer) or the DNA library of the
          pedigree index (germline); ``False`` otherwise
        * ``sex`` — ``"male"``, ``"female"``, or ``"unknown"``
        * ``donor_name`` — patient / family identifier
        * ``sample_name`` — bio-sample name
        * ``tissue_type`` — ``"tumor"``, ``"normal"``, or ``"unknown"``
        * ``is_background`` — whether the library is from a background data set; these are
          only included with ``include_background``
        * ``data_set`` — the name of the library's data set in the config
        * ``secondary_id`` — the secondary ID of the donor in the sample sheet

        Further sample sheet columns (``extra_infos`` of the library and its test sample,
        bio sample and bio entity) are added in snake_case, e.g. ``libraryKit`` becomes
        ``library_kit``; see :func:`_library_extra_columns`.
    """
    import pandas as pd
    from biomedsheets.shortcuts import GenericSampleSheet
    from biomedsheets.shortcuts.cancer import CancerCaseSheet, CancerCaseSheetOptions

    rows: list[dict] = []

    for info, _raw_sheet in zip(data_set_infos, sheets):
        if info.is_background and not include_background:
            continue

        extra_columns = _library_extra_columns(_raw_sheet) if _raw_sheet is not None else {}

        sheet_type = _sheet_type(info)

        # ── matched_cancer ───────────────────────────────────────────────────
        if sheet_type == "matched_cancer":
            try:
                csheet = CancerCaseSheet(
                    _raw_sheet,
                    options=CancerCaseSheetOptions(
                        allow_missing_normal=True, allow_missing_tumor=True
                    ),
                )
            except Exception as exc:
                print(
                    f"WARNING: could not build cancer library DataFrame: {exc}",
                    file=sys.stderr,
                )
                continue

            # Collect library names that belong to the primary tumor per donor
            # (first tumor in primary_sample_pairs).
            primary_tumor_lib_names: set[str] = set()
            for pair in csheet.primary_sample_pairs:
                if pair and pair.tumor_sample and pair.tumor_sample.dna_ngs_library:
                    primary_tumor_lib_names.add(pair.tumor_sample.dna_ngs_library.name)
                if pair and pair.tumor_sample and pair.tumor_sample.rna_ngs_library:
                    primary_tumor_lib_names.add(pair.tumor_sample.rna_ngs_library.name)

            for donor in csheet.donors:
                donor_name = donor.name
                for bio_sample in donor.bio_samples.values():
                    is_tumor = bio_sample.extra_infos.get("isTumor", False)
                    role = "tumor" if is_tumor else "normal"
                    sample_name = bio_sample.name
                    for ts in bio_sample.test_samples.values():
                        ext = ts.extra_infos.get("extractionType", "unknown").lower()
                        for lib in ts.ngs_libraries.values():
                            is_primary = lib.name in primary_tumor_lib_names if is_tumor else True
                            rows.append(
                                {
                                    "library_name": lib.name,
                                    "extraction_type": ext,
                                    "kind": "cancer",
                                    "role": role,
                                    "is_primary": is_primary,
                                    "sex": "unknown",
                                    "donor_name": donor_name,
                                    "sample_name": sample_name,
                                    "related_group": _group_prefix(info) + donor_name,
                                    "father_name": "0",
                                    "mother_name": "0",
                                    "disease_state": 0,
                                    "tissue_type": "tumor" if is_tumor else "normal",
                                    "is_background": info.is_background,
                                    "data_set": info.name,
                                    "secondary_id": donor.bio_entity.secondary_id,
                                    **extra_columns.get(lib.name, {}),
                                }
                            )

        # ── germline_variants (and generic) ──────────────────────────────────
        else:
            # All libraries, with the pedigrees of the germline view (joined by the data set's
            # ``pedigree_field``, if set)
            index_donor_names: set[str] = set()
            index_lib_names: set[str] = set()
            father_lib_names: set[str] = set()
            mother_lib_names: set[str] = set()
            lib_sex: dict[str, str] = {}
            donor_family: dict[str, str] = {}
            lib_father: dict[str, str] = {}
            lib_mother: dict[str, str] = {}
            lib_disease: dict[str, int] = {}

            all_ngs_libs = GenericSampleSheet(_raw_sheet).all_ngs_libraries
            germline_sheet = _germline_case_sheet(info, _raw_sheet)
            pedigrees = germline_sheet.cohort.pedigrees

            if pedigrees:
                # First pass: map donor PKs to their DNA library names
                pk_to_lib = {}
                for ped in pedigrees:
                    for d in ped.donors:
                        lib = getattr(d, "dna_ngs_library", None)
                        if lib:
                            pk_to_lib[getattr(d, "pk", None)] = lib.name

                families: dict[str, str] = {}
                for ped in pedigrees:
                    family = _family_name(info, ped)
                    if family in families:
                        raise ConfigError(
                            f"the pedigrees of {families[family]} and {_donor_name(ped.index)} in "
                            f"data set {info.name!r} are both named {family}; group them with "
                            f"pedigree_field: {FAMILY_FIELD} on the data set, or fix the sheet"
                        )
                    families[family] = _donor_name(ped.index)
                    # Index / proband: all libraries of the donor, so the RNA library too
                    if ped.index and ped.index.dna_ngs_library:
                        index_donor_names.add(_donor_name(ped.index))
                        index_lib_names.add(ped.index.dna_ngs_library.name)

                    # Build sets of PKs referenced as father / mother by children
                    pk_is_father: set = {
                        int(d.father_pk)
                        for d in ped.donors
                        if getattr(d, "father_pk", None) is not None
                    }
                    pk_is_mother: set = {
                        int(d.mother_pk)
                        for d in ped.donors
                        if getattr(d, "mother_pk", None) is not None
                    }

                    for d in ped.donors:
                        donor_family[_donor_name(d)] = family
                        lib = getattr(d, "dna_ngs_library", None)
                        if lib is None:
                            continue
                        pk = getattr(d, "pk", None)
                        sex_raw = (
                            d.extra_infos.get("sex", "unknown")
                            if hasattr(d, "extra_infos")
                            else "unknown"
                        )
                        lib_sex[lib.name] = sex_raw if sex_raw in ("male", "female") else "unknown"
                        if pk in pk_is_father:
                            father_lib_names.add(lib.name)
                        elif pk in pk_is_mother:
                            mother_lib_names.add(lib.name)

                        father_pk = getattr(d, "father_pk", None)
                        mother_pk = getattr(d, "mother_pk", None)
                        lib_father[lib.name] = pk_to_lib.get(
                            int(father_pk) if father_pk is not None else None, "0"
                        )
                        lib_mother[lib.name] = pk_to_lib.get(
                            int(mother_pk) if mother_pk is not None else None, "0"
                        )

                        is_affected = (
                            d.extra_infos.get("isAffected") if hasattr(d, "extra_infos") else None
                        )
                        if is_affected is True or str(is_affected).lower() == "affected":
                            lib_disease[lib.name] = 2
                        elif is_affected is False or str(is_affected).lower() == "unaffected":
                            lib_disease[lib.name] = 1
                        else:
                            lib_disease[lib.name] = 0

            for lib in all_ngs_libs:
                ts = lib.test_sample
                bs = getattr(ts, "bio_sample", None)
                bio_entity = getattr(bs, "bio_entity", None)

                # Handle both GenericSampleSheet (TestSample has extra_infos)
                # and GermlineCaseSheet (TestSampleShortcut wraps TestSample).
                ts_extra = getattr(ts, "extra_infos", None)
                if ts_extra is None and hasattr(ts, "test_sample"):
                    ts_extra = getattr(ts.test_sample, "extra_infos", {})
                else:
                    ts_extra = ts_extra or {}

                be_extra = getattr(bio_entity, "extra_infos", {}) if bio_entity else {}

                ext = ts_extra.get("extractionType", "unknown").lower()
                lib_name = lib.name
                donor_name = getattr(bio_entity, "name", "")
                sample_name = getattr(bs, "name", "")
                is_affected = be_extra.get("isAffected") or ts_extra.get("isAffected", None)
                # TSV sheets hold the words, JSON sheets booleans
                status = str(is_affected).lower() if is_affected is not None else "unknown"
                sex = lib_sex.get(lib_name, "unknown")

                # Only the index's DNA library is primary: the library that stands for the pedigree
                is_primary = lib_name in index_lib_names
                if donor_name in index_donor_names:
                    role = "index"
                elif lib_name in father_lib_names:
                    role = "father"
                elif lib_name in mother_lib_names:
                    role = "mother"
                elif is_affected is True or status == "affected":
                    role = "affected"
                elif is_affected is False or status == "unaffected":
                    role = "unaffected"
                elif not pedigrees:
                    role = "index"  # no pedigree information: each library stands alone
                    is_primary = True
                else:
                    role = "unknown"

                rows.append(
                    {
                        "library_name": lib_name,
                        "extraction_type": ext,
                        "kind": "germline",
                        "role": role,
                        "is_primary": is_primary,
                        "sex": sex,
                        "donor_name": donor_name,
                        "sample_name": sample_name,
                        "related_group": donor_family.get(
                            donor_name, _group_prefix(info) + donor_name
                        ),
                        "father_name": lib_father.get(lib_name, "0"),
                        "mother_name": lib_mother.get(lib_name, "0"),
                        "disease_state": lib_disease.get(lib_name, 0),
                        "tissue_type": "unknown",
                        "is_background": info.is_background,
                        "data_set": info.name,
                        "secondary_id": getattr(bio_entity, "secondary_id", ""),
                        **extra_columns.get(lib_name, {}),
                    }
                )

    extra_column_names = sorted({column for row in rows for column in row} - set(_LIBRARY_COLUMNS))
    columns = _LIBRARY_COLUMNS + extra_column_names
    df = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)
    if relationships:
        df = resolve_relationships(df, relationships)
    return df


#: The columns of the donor table of :func:`build_donor_dataframe`
_DONOR_COLUMNS = [
    "donor_name",
    "father_name",
    "mother_name",
    "sex",
    "library_name",
    "pedigree",
    "data_set",
    "is_background",
]


def _donor_name(donor) -> str:
    """The name of a donor as the library table names it: that of its bio entity."""
    return getattr(getattr(donor, "bio_entity", None), "name", None) or donor.name


# ---------------------------------------------------------------------------
# Related groups: their names, entities and libraries
# ---------------------------------------------------------------------------


#: The custom field that names a family on germline sheets, such as those of ``snappy pull-sheet``
FAMILY_FIELD = "familyId"

#: The prefix of a related group that is named after a donor, by sheet type
_GROUP_PREFIXES = {"matched_cancer": "PAT_", "germline_variants": "FAM_", "generic": ""}


def _sheet_type(info) -> str:
    return getattr(info, "sheet_type", "") or str(getattr(info, "type", ""))


def _germline_case_sheet(info, raw_sheet):
    """The biomedsheets germline view of a sheet, with the pedigrees joined by the data set's ``pedigree_field``.
    Raises ``ConfigError`` when a donor lacks the field."""
    from biomedsheets.shortcuts import GermlineCaseSheet
    from biomedsheets.shortcuts.germline import UndefinedFieldException

    try:
        return GermlineCaseSheet(raw_sheet, **(getattr(info, "pedigree_field_kwargs", None) or {}))
    except UndefinedFieldException as exc:
        raise ConfigError(f"pedigree_field of data set {info.name!r}: {exc}") from exc


def _group_prefix(info) -> str:
    """The prefix of the related groups of a data set that are named after a donor."""
    return _GROUP_PREFIXES.get(_sheet_type(info), "")


def _family_name(info, pedigree) -> str:
    """Return the name of a pedigree, the ``related_group`` of its members.
    It is the value of the data set's ``pedigree_field``, else of the ``familyId`` custom field, over the donors that have one.
    Without any, it is the index donor's name with the prefix of the sheet type: ``FAM_`` on germline sheets, none on generic ones.
    Raises ``ConfigError`` when the donors disagree: the pedigree then comes from the parent columns, and the field names another family."""
    field = (getattr(info, "pedigree_field_kwargs", None) or {}).get("join_by_field")
    name_field = field or FAMILY_FIELD
    values = {
        str(value)
        for donor in pedigree.donors
        if (value := getattr(donor, "extra_infos", {}).get(name_field)) not in (None, "", ".")
    }
    index = _donor_name(pedigree.index)
    if len(values) > 1:
        raise ConfigError(
            f"the pedigree of {index} has several values of {name_field}: "
            f"{', '.join(sorted(values))}; group the pedigrees by the field with pedigree_field: "
            f"{name_field} on the data set, or fix the sheet"
        )
    if field and not values:
        raise ConfigError(
            f"{field} of data set {info.name!r} has no value for the pedigree of {index}"
        )
    return values.pop() if values else _group_prefix(info) + index


def entity_libraries(df, group_by: str | None) -> dict[str, str]:
    """Return the output entities of the libraries ``df`` with their primary library.
    With ``group_by == "related"`` each ``related_group`` is an entity, named like the group, and its primary library is the ``is_primary`` row, else the first one; otherwise every library is an entity of its own."""
    if df.empty:
        return {}
    if group_by == "related":
        entities = {}
        for name, group in df.groupby("related_group"):
            primary = group[group["is_primary"]]
            entities[name] = (primary if not primary.empty else group).iloc[0]["library_name"]
        return entities
    return {name: name for name in df["library_name"].unique()}


def _check_groups(df) -> None:
    """Raise ``ConfigError`` for a related group whose libraries ``df`` span several data sets,
    since the group's files are named after it."""
    if df.empty or not {"related_group", "data_set"} <= set(df.columns):
        return
    data_sets = df.groupby("related_group")["data_set"].agg(lambda names: sorted(set(names)))
    if colliding := {name: sets for name, sets in data_sets.items() if len(sets) > 1}:
        names = ", ".join(f"{name} ({', '.join(sets)})" for name, sets in colliding.items())
        first_data_set = next(iter(colliding.values()))[0]
        raise ConfigError(
            f"related groups of the same name in several data sets: {names}. A task names the "
            "files of a related group after it, so it can only take one group of a name. A family "
            "sequenced in batches belongs into one sheet: one data set, whose search_paths list "
            "the folders of all batches. For unrelated groups, select the libraries of one data "
            f'set with library_selection, for example library_selection: "data_set == '
            f"'{first_data_set}'\", or use one task per data set."
        )


def build_donor_dataframe(data_set_infos, sheets, include_background=False):
    """Build a DataFrame with one row per donor of the germline sample sheets, sequenced or
    not: the parents as donor names (``"0"`` for none), ``sex``, the donor's primary DNA
    library (``""`` without one) and ``pedigree``, the name of its pedigree (the ``related_group``
    of its members). Cancer sheets have no rows: their donors are all sequenced."""
    import pandas as pd

    rows: list[dict] = []
    for info, raw_sheet in zip(data_set_infos, sheets):
        if info.is_background and not include_background:
            continue
        sheet_type = _sheet_type(info)
        if sheet_type == "matched_cancer" or raw_sheet is None:
            continue
        germline_sheet = _germline_case_sheet(info, raw_sheet)
        for pedigree in germline_sheet.cohort.pedigrees:
            family = _family_name(info, pedigree)
            for donor in pedigree.donors:
                parents = []
                for pk in (donor.father_pk, donor.mother_pk):
                    parent = pedigree.pk_to_donor.get(str(pk)) if pk is not None else None
                    parents.append(_donor_name(parent) if parent else "0")
                sex = (
                    donor.extra_infos.get("sex", "unknown")
                    if hasattr(donor, "extra_infos")
                    else "unknown"
                )
                rows.append(
                    {
                        "donor_name": _donor_name(donor),
                        "father_name": parents[0],
                        "mother_name": parents[1],
                        "sex": sex if sex in ("male", "female") else "unknown",
                        "library_name": donor.dna_ngs_library.name if donor.dna_ngs_library else "",
                        "pedigree": family,
                        "data_set": info.name,
                        "is_background": info.is_background,
                    }
                )
    return pd.DataFrame(rows, columns=_DONOR_COLUMNS)


def resolve_relationships(df, relationships: dict):
    """Add relationship columns to the library DataFrame.

    For each relationship definition, a new column is added whose value is
    looked up from related rows sharing the same ``via`` column value and
    matching the ``target`` query expression.

    Parameters
    ----------
    df:
        Library DataFrame produced by :func:`build_library_dataframe`.
    relationships:
        Dict mapping relationship names to :class:`RelationshipDefinition`
        instances (or dicts with ``via``, ``target``, ``column``, ``many``
        keys).

    Returns
    -------
    pandas.DataFrame
        Copy of *df* with additional relationship columns.
    """
    if df.empty or not relationships:
        return df

    result = df.copy()

    for rel_name, rel_def in relationships.items():
        # Accept both model instances and plain dicts.
        via = rel_def.via if hasattr(rel_def, "via") else rel_def["via"]
        target = rel_def.target if hasattr(rel_def, "target") else rel_def["target"]
        column = (
            rel_def.column if hasattr(rel_def, "column") else rel_def.get("column")
        ) or rel_name
        many = rel_def.many if hasattr(rel_def, "many") else rel_def.get("many", False)

        via_columns = [via] if isinstance(via, str) else list(via)
        for via_column in via_columns:
            if via_column not in result.columns:
                raise ConfigError(
                    f"Relationship {rel_name!r}: join column {via_column!r} not found in "
                    f"DataFrame. Available columns: {list(result.columns)}"
                )
        # One join key per row; several ``via`` columns join on all of them
        via_keys = result[via_columns].astype(str).agg("/".join, axis=1)

        # Build a lookup: for each unique ``via`` value, find matching rows.
        # We query the full DataFrame for each group so that the target
        # expression can reference any column.
        matched: dict[str, list[str]] = {}
        for via_val, group in result.groupby(via_keys):
            try:
                hits = group.query(target)
            except Exception as exc:
                raise ConfigError(
                    f"Relationship {rel_name!r}: target query {target!r} failed: {exc}"
                ) from exc
            matched[str(via_val)] = hits["library_name"].tolist()

        if many:
            result[column] = via_keys.map(lambda v: matched.get(str(v), []))
        else:

            def _single_match(via_value):
                matches = matched.get(str(via_value), [])
                if len(matches) > 1:
                    raise ConfigError(
                        f"Relationship {rel_name!r}: expected single match for "
                        f"via value {via_value!r}, got {len(matches)}: {matches}. "
                        f"Use 'many: true' if multiple matches are expected."
                    )
                return matches[0] if matches else ""

            result[column] = via_keys.map(_single_match)

    return result


def select_libraries(
    df, selection: str | None, default: str | None = None, requirement: str | None = None
):
    """Return the libraries of ``df`` that a task takes.

    These are the rows of the step's ``requirement`` (what it can process at all) that are
    also rows of ``selection`` (the task's ``library_selection``), else those of the
    step's ``default``, else those of the default of each row's sheet kind
    (``_LIBRARY_SELECTION_DEFAULTS``); an empty expression takes all rows. Raises ``ConfigError``
    for an invalid expression and for libraries of the same name in several data sets.
    """
    import pandas as pd

    if df.empty:
        return df
    if requirement:
        df = _query_libraries(df, requirement)
    expr = selection if selection is not None else default
    if expr is None:
        df = pd.concat(
            [
                subset if subset.empty else _query_libraries(subset, kind_default)
                for kind, kind_default in _LIBRARY_SELECTION_DEFAULTS.items()
                for subset in [df[df["kind"] == kind]]
            ]
        )
    elif expr:
        df = _query_libraries(df, expr)
    if colliding := _colliding_libraries(df):
        raise ConfigError(colliding_libraries_message(colliding))
    _check_groups(df)
    return df


def _query_libraries(df, expr: str):
    try:
        return df.query(expr)
    except Exception as exc:
        raise ConfigError(f"library_selection expression {expr!r} is invalid: {exc}") from exc


def _colliding_libraries(df) -> dict[str, list[str]]:
    """Return library name -> data sets, for the names of several libraries in ``df``."""
    duplicated = df[df["library_name"].duplicated(keep=False)]
    return {name: list(group["data_set"]) for name, group in duplicated.groupby("library_name")}


def colliding_libraries_message(colliding: dict[str, list[str]]) -> str:
    """Return the error for libraries of the same name from several data sets, with the fix."""
    names = ", ".join(f"{name} ({', '.join(data_sets)})" for name, data_sets in colliding.items())
    first_data_set = next(iter(colliding.values()))[0]
    return (
        f"libraries of the same name in several data sets: {names}. A task names its files by "
        "library, so it can only take one library of a name. Select the libraries of one data "
        f'set with library_selection, for example library_selection: "data_set == '
        f"'{first_data_set}'\", or use one task per data set."
    )
