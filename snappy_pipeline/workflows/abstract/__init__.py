# -*- coding: utf-8 -*-
"""Base classes for the actual pipeline steps"""

import contextlib
import datetime
import logging
import os
import os.path
import sys
import tempfile
import typing
from collections import OrderedDict
from collections.abc import MutableMapping
from dataclasses import dataclass
from fnmatch import fnmatch
from functools import lru_cache
from io import StringIO
from typing import Any, Callable

import pydantic
import ruamel.yaml as ruamel_yaml
from biomedsheets import io_tsv
from biomedsheets.io import SheetBuilder, json_loads_ordered
from biomedsheets.models import SecondaryIDNotFoundException
from biomedsheets.naming import NAMING_SCHEMES, name_generator_for_scheme
from biomedsheets.ref_resolver import RefResolver
from biomedsheets.shortcuts import (
    ShortcutSampleSheet,
    donor_has_dna_ngs_library,
    write_pedigree_to_ped,
    write_pedigrees_to_ped,
)
from snakemake.api import Workflow
from snakemake.io import touch
from snakemake.iocontainers import InputFiles, OutputFiles, Wildcards

from snappy_pipeline.base import (
    MissingConfiguration,
    UnsupportedActionException,
    merge_kwargs,
    print_config,
    print_sample_sheets,
)
from snappy_pipeline.find_file import FileSystemCrawler, PatternSet
from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract.pedigree import append_pedigree_to_ped
from snappy_pipeline.workflows.abstract.protocol import DataSignature, ExpectedPathSchema
from snappy_wrappers.resource_usage import ResourceUsage

#: String constant with bash command for redirecting stderr to ``{log}`` file
STDERR_TO_LOG_FILE = r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file and enable printing executed commands
exec 2> >(tee -a "{log}")
set -x
# -----------------------------------------------------------------------------

""".lstrip()


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


class ImplementationUnavailableError(NotImplementedError):
    """Raised when a function that is to be overridden optionally is called

    This is provided as an alternative to ``NotImplementedError`` as the Python linters warn if
    a class does not override functions throwing ``NotImplementedError``.
    """


Inputs: typing.TypeAlias = InputFiles | dict[str, Any]
Outputs: typing.TypeAlias = OutputFiles | dict[str, Any]
Threads: typing.TypeAlias = int | None
Attempt: typing.TypeAlias = int | None


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
        """Return the resource usage for the given action."""
        if action not in self.actions:
            raise ValueError(f"Invalid {action} not in {self.actions}")
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
            else:
                return getattr(resource_usage, resource_name)

        return _get_resource

    def get_args(self, action: str) -> Inputs | Callable[[Wildcards], Inputs]:
        """Return args for the given action of the sub step"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_input_files(self, action: str) -> Inputs | Callable[[Wildcards], Inputs]:
        """Return input files for the given action of the sub step"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_output_files(self, action: str) -> Outputs:
        """Return output files for the given action of the sub step and"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_log_file(self, action: str) -> Outputs:
        """Return path to log file

        The default implementation tries to call ``self._get_log_files()`` and in the case of
        this function returning a ``dict``, augments it with paths to MD5 files.
        """
        if hasattr(self, "_get_log_file"):
            inner = self._get_log_file(action)
            if not isinstance(inner, dict):
                return inner
            else:
                result = {}
                for k, v in inner.items():
                    result[k] = v
                    if not k.endswith("_md5") and (k + "_md5") not in inner:
                        result[k + "_md5"] = v + ".md5"
                return result
        else:
            raise NotImplementedError(
                "Log file name generation not implemented!"
            )  # pragma: no cover

    def get_shell_cmd(self, action: str, wildcards: Wildcards) -> str:  # NOSONAR
        """Return shell command for the given action of the sub step and the given wildcards"""
        raise ImplementationUnavailableError(
            "Override this method before calling it!"
        )  # pragma: no cover

    def run_locally(self, action: str, wildcards: Wildcards) -> str:  # NOSONAR
        """Runs a function locally for the given action of the sub step and the given wildcards"""
        raise ImplementationUnavailableError(
            "Override this method before calling it!"
        )  # pragma: no cover

    def run(self, action: str, wildcards: Wildcards):  # NOSONAR
        """Run the sub steps action action's code with the given wildcards"""
        raise ImplementationUnavailableError(
            "Override this method before calling it!"
        )  # pragma: no cover

    def check_config(self):
        """Check configuration, raise ``ConfigurationMissing`` on problems

        Override in sub classes.

        :raises:MissingConfiguration: on missing configuration
        """


class WritePedigreeStepPart(BaseStepPart):
    """Write out pedigree file for primary DNA sample given the index NGS library name"""

    #: Step name
    name = "write_pedigree"

    #: Class available actions
    actions = ("run",)

    def __init__[P: BaseStep](
        self, parent: P, require_dna_ngs_library: bool = False, only_trios: bool = False
    ):
        super().__init__(parent)
        #: Whether to prevent writing out of samples with out NGS library.
        self.require_dna_ngs_library = require_dna_ngs_library
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            if require_dna_ngs_library:
                for name, pedigree in sheet.index_ngs_library_to_pedigree.items():
                    pedigree = pedigree.with_filtered_donors(donor_has_dna_ngs_library)
                    if only_trios:
                        in_trio = set()
                        for donor in pedigree.donors:
                            if donor.father and donor.mother:
                                in_trio |= {
                                    donor.name,
                                    donor.father.name,
                                    donor.mother.name,
                                }
                        if not any((donor.name in in_trio for donor in pedigree.donors)):
                            continue  # ignore empty pedigree post filtration
                        pedigree = pedigree.with_filtered_donors(
                            lambda donor: donor.name in in_trio
                        )
                    self.index_ngs_library_to_pedigree[name] = pedigree
            else:
                self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Returns function returning input files.

        Returns a dict with entry ``"bam"`` mapping to list of input BAM files.  This list will
        be empty if the parent step does not define an ``"ngs_mapping"`` workflow.
        """
        self._validate_action(action=action)

        @listify
        def get_input_files(wildcards):
            if "ngs_mapping" not in self.parent.modules:
                return  # early exit
            # Get shortcut to NGS mapping sub workflow
            ngs_mapping = self.parent.modules["ngs_mapping"]
            # Get names of primary libraries of the selected pedigree.  The pedigree is selected
            # by the primary DNA NGS library of the index.
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            if not pedigree.index or not pedigree.index.dna_ngs_library:
                msg = "INFO: pedigree without index (names: {})"  # pragma: no cover
                donor_names = list(sorted(d.name for d in pedigree.donors))
                print(msg.format(donor_names), file=sys.stderr)  # pragma: no cover
                return
            tpl = "output/{library_name}/out/{library_name}{ext}"
            for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
                library_name = donor.dna_ngs_library.name
                path = tpl.format(
                    library_name=library_name,
                    ext=".bam",
                    **wildcards,
                )
                yield ngs_mapping(path)

        return get_input_files

    def get_output_files(self, action):
        self._validate_action(action=action)
        return "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

    # @listify
    def get_result_files(self):
        # tpl = self.get_output_files("run")
        # for sheet in getattr(self.parent, "shortcut_sheets", []):
        #     for index_ngs_library in sheet.index_ngs_library_to_pedigree.keys():
        #         yield tpl.format(index_ngs_library=index_ngs_library)
        return []

    def run(self, wildcards: Wildcards, output: OutputFiles):
        """Write out the pedigree information

        :param wildcards: Snakemake wildcards associated with rule (unused).
        :type wildcards: snakemake.io.Wildcards

        :param output: Snakemake output associated with rule.
        :type output: snakemake.io.Namedlist
        """
        if wildcards.index_ngs_library == "whole_cohort":
            write_pedigrees_to_ped(self.index_ngs_library_to_pedigree.values(), str(output))
        else:
            write_pedigree_to_ped(
                self.index_ngs_library_to_pedigree[wildcards.index_ngs_library],
                str(output),
            )


class WritePedigreeSampleNameStepPart(WritePedigreeStepPart):
    """
    Class contains method to write pedigree file for primary DNA sample given the index
    NGS library name.It will create pedigree information based sole on sample name,
    example 'P001' instead of 'P001-N1-DNA1-WGS1'.
    """

    #: Step name
    name = "write_pedigree_with_sample_name"

    def __init__(self, *args, **kwargs):
        WritePedigreeStepPart.__init__(self, *args, **kwargs)

    def run(self, wildcards, output):
        """Write out the pedigree information

        :param wildcards: Snakemake wildcards associated with rule (unused).
        :type wildcards: snakemake.io.Wildcards

        :param output: Snakemake output associated with rule.
        :type output: snakemake.io.Namedlist
        """
        append_pedigree_to_ped(
            pedigree=self.index_ngs_library_to_pedigree[wildcards.index_ngs_library],
            output_path=str(output),
        )


class LinkOutStepPart(BaseStepPart):
    """Generically link out

    This is for output files that are created unconditionally, i.e., for output files where the
    output name is the same as for the work file.
    """

    name = "link_out"

    def __init__(self, parent, disable_patterns=None):
        super().__init__(parent)
        self.base_pattern_out = "output/{path}/{file}.{ext}"
        self.base_path_out = self.base_pattern_out.replace(",[^/]", "")
        self.base_path_in = "work/{path}/{file}.{ext}"
        #: Patterns for disabling linking out to.  This is useful/required when there is a
        #: specialized link out step part, e.g., for the case of alignment where realignment is
        #: performed or not, depending on the configuration.
        self.disable_patterns = list(disable_patterns or [])

    def get_input_files(self, action):
        """Return input file pattern"""

        if not self.disable_patterns:
            return self.base_path_in

        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""

        def input_function(wildcards):
            """Helper wrapper function"""
            result = self.base_path_in.format(**wildcards)
            for pattern in self.disable_patterns:
                if fnmatch(result, pattern):
                    raise ValueError("Blocking match...")
            return task_prefix + result

        return input_function

    def get_output_files(self, action):
        """Return output file pattern"""
        assert action == "run", "Unsupported action"
        return self.base_pattern_out

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out"""
        assert action == "run", "Unsupported action"
        return "test -h {output[0]} || ln -sr {input[0]} {output[0]}"

    def run_locally(self, action, wildcards):
        assert action == "run", "Unsupported action"
        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""
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


@lru_cache()
def _cached_yaml_round_trip_load_str(str_value):
    """Cached reading of YAML ``str`` objects."""
    yaml = ruamel_yaml.YAML()
    return yaml.load(StringIO(str_value))


class DataSetInfo:
    """Information on a DataSet"""

    def __init__(
        self,
        name,
        sheet_path,
        base_paths,
        search_paths,
        search_patterns,
        sheet_type,
        is_background,
        naming_scheme,
        mixed_se_pe,
        sodar_uuid,
        sodar_title,
        pedigree_field=None,
    ):
        """Constructor.

        :param name: Name of the data set.
        :type name: str

        :param sheet_path: Path to sheet file that should be loaded.
        :type sheet_path: str

        :param base_paths:  All base paths of all configuration, to look for ``sheet_path``.
        :type base_paths: list

        :param search_paths: Search paths for the files in the sample sheet.
        :type search_paths: list

        :param search_patterns: Search patterns. Example: "{left: '**/*_R1_*.fastq.gz',
        right: '**/*_R2_*.fastq.gz'}".
        :type search_patterns: List[Dict[str, str]]

        :param sheet_type: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will
        attempt to load from sheet.
        :type sheet_type: str

        :param is_background: Whether the data set info is to be used only for background.
        :type is_background: bool

        :param naming_scheme: Selected naming schema: either 'secondary_id_pk' or
        'only_secondary_id'.
        :type naming_scheme: str

        :param mixed_se_pe: Whether mixing SE and PE data sets is allowed.

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
        #: All base paths of all configuration, to look for ``sheet_path``
        self.base_paths = base_paths
        #: Search paths for the files in the sample sheet
        self.search_paths = list(search_paths)
        #: Search patterns
        self.search_patterns = search_patterns
        #: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will attempt to load
        # from sheet.
        self.sheet_type = sheet_type
        #: Whether the data set info is to be used only for background
        self.is_background = is_background
        #: Selected naming schema
        if naming_scheme not in NAMING_SCHEMES:
            raise ValueError("Invalid naming scheme: {}".format(naming_scheme))  # pragma: no cover
        self.naming_scheme = naming_scheme
        #: Whether mixing SE and PE data sets is allowed.
        self.mixed_se_pe = mixed_se_pe
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
        for base in self.base_paths:
            fname = os.path.join(base, self.sheet_path)
            if not os.path.exists(fname):
                continue
            if self.sheet_path.endswith(".tsv"):
                if self.sheet_type == "matched_cancer":
                    return self._set_is_background(
                        _cached_read_cancer_tsv_sheet(fname, self.sheet_path, self.naming_scheme),
                        self.is_background,
                    )
                elif self.sheet_type == "germline_variants":
                    return self._set_is_background(
                        _cached_read_germline_tsv_sheet(fname, self.sheet_path, self.naming_scheme),
                        self.is_background,
                    )
                elif self.sheet_type == "generic":
                    return self._set_is_background(
                        _cached_read_generic_tsv_sheet(fname, self.sheet_path, self.naming_scheme),
                        self.is_background,
                    )
                else:
                    raise ValueError("Invalid sheet type {}".format(self.sheet_type))
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
        # Raise problem if we could not find the sample sheet
        raise ValueError(  # pragma: no cover
            "Could not find sample sheet file {} in the lookup paths {}".format(
                self.sheet_path, self.base_paths
            )
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
class DataSearchInfo:
    """Data search information - simplified version of ``DataSetInfo``."""

    sheet_path: str
    base_paths: list
    search_paths: list
    search_patterns: list
    mixed_se_pe: bool


class BaseStep:
    """Base class for the pipeline steps

    Each pipeline step is a Snakemake workflow
    """

    #: Override with step name
    name: str

    #: DataSignatures this workflow step consumes
    consumes: dict[DataSignature, bool] = {}

    #: DataSignatures this workflow step produces
    produces: list[DataSignature] = []

    def __init_subclass__(cls, **kwargs):
        """Warn when a ``BaseStep`` subclass does not declare a ``produces`` contract."""
        super().__init_subclass__(**kwargs)
        # Only warn for concrete leaf classes that inherit produces=[] from BaseStep
        own_produces = cls.__dict__.get("produces")
        if own_produces is None and not cls.produces:
            import warnings as _warnings

            _warnings.warn(
                f"Transition Guidance: '{cls.__name__}' should declare a 'produces' list "
                "for full contract compliance.",
                DeprecationWarning,
                stacklevel=2,
            )

    #: Override with the sheet shortcut class to use
    sheet_shortcut_class: type[ShortcutSampleSheet]

    #: Override with arguments to pass into sheet shortcut class constructor
    sheet_shortcut_args = None

    #: Override with keyword arguments to pass into sheet shortcut class
    #: constructor
    sheet_shortcut_kwargs = None

    @classmethod
    def default_config_yaml(cls):
        """Override this function for providing default configuration

        The configuration should be a YAML fragment. Your configuration should define a top-level
        key starting with '_' and then consist of the name of the schema, e.g.,
        '_ngs_mapping_schema'. Your default configuration is then merged into the main
        configuration where the main configuration takes precedence.

        Example: ::

            def default_config_yaml(self):
                return textwrap.dedent(""\"
                    schema_config:
                      ngs_mapping:
                        max_threads: 16
                ""\").lstrip()))

        Return ``None`` for no default configuration.

        You can also return an iterable of configurations, these will be merged in the order given
        (earlier ones will be overwritten by later ones).  This is useful if your schema needs
        configuration for a later one.
        """
        return ""  # pragma: no cover

    #: Override with the Pydantic model class for configuration validation
    config_model_class: type[SnappyStepModel]

    @classmethod
    def get_output_paths(
        cls, signature: "DataSignature | None" = None, **kwargs
    ) -> "dict[str, str]":
        """Return **local** output paths for *signature*.

        Override in concrete workflow classes to expose this step's outputs to downstream
        consumers via :meth:`get_upstream_paths`.

        Arguments:
            signature: The :class:`~snappy_pipeline.workflows.abstract.protocol.DataSignature`
                requested by the consumer.  Implementations should check
                ``signature.satisfies(...)`` and raise ``ValueError`` for unsupported
                signatures.
            **kwargs: Caller-supplied identifiers (e.g. ``library_name``, ``sample_name``)
                used to parametrise the returned path strings.  Implementations may use
                ``{library_name}``-style format strings when the identifiers are omitted.

        Returns:
            A ``dict[str, str]`` mapping logical output key (e.g. ``"bam"``, ``"vcf"``) to a
            **local** path string (relative to the step's own ``output/`` directory, e.g.
            ``"output/{library_name}/out/{library_name}.bam"``).

        Raises:
            NotImplementedError: When the concrete subclass has not overridden this method.
            ValueError: When *signature* is not supported by this workflow.
        """
        raise NotImplementedError(
            f"'{cls.__name__}' must implement 'get_output_paths' to act as an upstream provider."
        )

    def __init__(
        self,
        workflow: Workflow,
        config: MutableMapping[str, Any],
        config_lookup_paths: tuple[str, ...],
        config_paths: tuple[str, ...],
        work_dir: str,
        *,
        task_name: str | None = None,
        previous_steps: tuple[type[typing.Self], ...] | None = None,
    ):
        self.step_name = self.__class__.name
        self.config_paths = config_paths
        self.work_dir = work_dir
        self.previous_steps = tuple(previous_steps or [])
        self.workflow = workflow
        self.modules = {}

        try:
            from snappy_pipeline.workflow_model import ConfigModel

            self.w_config: ConfigModel = ConfigModel(**config)
        except pydantic.ValidationError as ve:
            raise ve

        # 1. Look up the task name injected by the orchestrator
        req_task_name = config.get("__task_name__") or task_name or self.step_name
        self.task = next((t for t in self.w_config.tasks if t.name == req_task_name), None)

        # 2. If the task doesn't match our step type (i.e. we are a submodule being
        # blindly initialized by a parent's boilerplate Snakefile), ignore it and
        # grab the actual config for our step type.
        if not self.task or self.task.step != self.step_name:
            self.task = next((t for t in self.w_config.tasks if t.step == self.step_name), None)

        if not self.task:
            raise ValueError(f"No task configuration found for step '{self.step_name}'.")

        self.task_name = self.task.name
        self.logger = logging.getLogger(self.task_name)

        print(f"\n[DEBUG] Initializing step '{self.step_name}' for task '{self.task_name}'")

        # Validate from mapping input explicitly to ensure nested coercion is applied consistently.
        self.config = self.config_model_class.model_validate(self.task.config)
        self.depends_on = getattr(self.config, "depends_on", None)

        self.config_lookup_paths = list(config_lookup_paths)
        self.sub_steps: dict[str, BaseStepPart] = {}
        self.data_set_infos = list(self._load_data_set_infos())
        self.sheets = [info.sheet for info in self.data_set_infos]

        sheet_kwargs_list = [
            merge_kwargs(
                first_kwargs=self.sheet_shortcut_kwargs,
                second_kwargs=info.pedigree_field_kwargs,
            )
            for info in self.data_set_infos
        ]

        self.shortcut_sheets = []
        klass = self.__class__.sheet_shortcut_class
        for sheet, kwargs in zip(self.sheets, sheet_kwargs_list):
            kwargs = kwargs or {}
            kwargs = {k: v for k, v in kwargs.items() if k in klass.supported_kwargs}
            self.shortcut_sheets.append(
                klass(sheet, *(self.__class__.sheet_shortcut_args or []), **kwargs)
            )

        self._setup_hooks()
        self._check_config()

        config_string = self.config.model_dump_yaml(by_alias=True)

        _config = _cached_yaml_round_trip_load_str(config_string)
        config.update(_config)

    def get_task_config(self, name: str) -> SnappyStepModel:
        """Retrieve the typed configuration model of an upstream task based on dependency resolution."""

        # If the requested name matches this instance's step type or task name, return own config.
        if name == self.name or name == getattr(self, "task_name", ""):
            return self.config

        # Resolve via config-model typed dependency mapping, then literal name as fallback.
        dep_target = getattr(self.depends_on, name, None) if self.depends_on is not None else None
        if dep_target is None and self.depends_on is not None:
            # Check if any resolved dependency task has a step that matches `name`
            for dep_field, dep_val in self.depends_on.model_dump().items():
                if isinstance(dep_val, str) and dep_val:
                    t = next((tk for tk in self.w_config.tasks if tk.name == dep_val), None)
                    if t and t.step == name:
                        dep_target = dep_val
                        break
        target_task_name = dep_target or name

        # Find the task in the global config
        task = next((t for t in self.w_config.tasks if t.name == target_task_name), None)

        if not task:
            matching_steps = [t for t in self.w_config.tasks if t.step == target_task_name]
            if len(matching_steps) == 1:
                task = matching_steps[0]
            elif len(matching_steps) > 1:
                raise ValueError(
                    f"Ambiguous dependency: '{target_task_name}' matches multiple tasks by step type. "
                    f"Please explicitly map it in the 'depends_on' block for task '{self.task_name}'."
                )

        if not task:
            raise ValueError(
                f"Task '{target_task_name}' (resolved from '{name}') not found in configuration."
            )

        # Instantiate and return its strictly typed config model
        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        wf_class = WORKFLOW_REGISTRY.get(task.step)
        if not wf_class:
            raise ValueError(
                f"Workflow class for step '{task.step}' not found in WORKFLOW_REGISTRY."
            )

        return wf_class.config_model_class.model_validate(task.config)

    def get_preprocessed_path(self) -> str:
        """Return a preprocessed FASTQ directory from configured RAW dependencies.

        Resolution order:

        1. The ``link_in`` dependency, if configured, using its explicit ``path`` value.
        2. Any other configured ``depends_on`` field annotated with ``DataSignature(DataType.RAW)``,
           interpreted as an in-pipeline task that exposes FASTQs under
           ``<upstream_task_name>/output`` (e.g. ``adapter_trimming``).

        Returns an empty string if no matching RAW provider dependency is configured.
        """
        if self.depends_on is None:
            return ""

        # Prefer explicit external path from link_in, if available.
        link_in_task_name = getattr(self.depends_on, "link_in", "")
        if link_in_task_name:
            try:
                upstream_config = self.get_task_config("link_in")
                explicit_path = getattr(upstream_config, "path", "") or ""
                if explicit_path:
                    return explicit_path
            except Exception:
                pass

        # Fallback: any configured RAW dependency task with standard output layout.
        for field_name, field_info in type(self.depends_on).model_fields.items():
            dep_task_name = getattr(self.depends_on, field_name, "")
            if not dep_task_name:
                continue
            is_raw_dep = any(
                isinstance(meta, DataSignature) and meta.type.value == "raw"
                for meta in field_info.metadata
            )
            if not is_raw_dep:
                continue
            if field_name == "link_in":
                continue
            return f"{dep_task_name}/output"

        return ""

    def _setup_hooks(self):
        """Setup Snakemake workflow hooks for start/end/error"""

        # In the following, the "log" parameter to the handler functions is set to "_" as we
        # don't use them
        def on_start(_):
            """Print configuration and sample sheets on start"""
            verbose = False
            if verbose:
                # Print configuration back to the user after merging workflow step-specific
                # configuration
                print_config(self.config, file=sys.stderr)
                # Print sample sheets to the user
                print_sample_sheets(self, file=sys.stderr)

        def on_error(_):
            """Error handler, print message"""
            msg = "Oh no! Something went wrong."
            print("\n" + "*" * len(msg), file=sys.stderr)
            print(msg, file=sys.stderr)
            print("*" * len(msg) + "\n", file=sys.stderr)

        def on_success(_):
            """Success handler, print message"""
            msg = "All done; have a nice day!"
            print("\n" + "*" * len(msg), file=sys.stderr)
            print(msg, file=sys.stderr)
            print("*" * len(msg) + "\n", file=sys.stderr)

        self.workflow.onstart(on_start)
        self.workflow.onerror(on_error)
        self.workflow.onsuccess(on_success)

    def _check_config(self):
        """Internal method, checks step and sub step configurations"""
        self.check_config()
        for step in self.sub_steps.values():
            step.check_config()

    def check_config(self):
        """Check ``self.w_config``, raise ``ConfigurationMissing`` on problems

        Override in sub classes.

        :raises:MissingConfiguration: on missing configuration
        """

    def ensure_w_config(self, config_keys, msg, e_class=MissingConfiguration):
        """Check parameters in configuration.

        Method ensures required configuration setting are present in the provided configuration;
        if not, it raises exception.

        :param config_keys: List of strings with all keys that must be present in the configuration
        for a given step of the analysis to be performed.
        :type config_keys: tuple

        :param msg: Message to be used in case of exception.
        :type msg: str

        :param e_class: Preferred exception class to be raised in case of error.
        Default: MissingConfiguration.
        :type e_class: class
        """
        # Initialise variables
        so_far = []
        handle = self.w_config

        # Check if configuration is empty
        if not handle:
            tpl = 'Empty configuration ("{full_path}"): {msg}'.format(
                full_path="/".join(config_keys), msg=msg
            )
            raise e_class(tpl)

        # Iterate over required configuration keys
        for entry in config_keys:
            # Check if keys are present in config dictionary
            if (handle := (getattr(handle, entry, None) or handle.get(entry, None))) is not None:
                so_far.append(entry)
            else:
                tpl = 'Missing configuration ("{full_path}", got up to "{so_far}"): {msg}'.format(
                    full_path="/".join(config_keys), so_far="/".join(so_far), msg=msg
                )
                raise e_class(tpl)

    def register_sub_step_classes(
        self, classes: tuple[type[BaseStepPart] | tuple[type[BaseStepPart], Any], ...]
    ):
        """Register an iterable of sub step classes

        Initializes objects in ``self.sub_steps`` dict
        """
        for pair_or_class in classes:
            try:
                klass, args = pair_or_class
            except TypeError:
                klass = pair_or_class
                args = ()
            obj = klass(self, *args)
            # obj.check_config()
            self.sub_steps[klass.name] = obj

    def get_upstream_dir(self, requirement: DataSignature) -> str:
        """
        Resolves the physical output directory for a required dependency
        based on the 'depends_on' configuration mapping.
        """
        upstream_task_name = (
            getattr(self.depends_on, requirement.type.value, None)
            if self.depends_on is not None
            else None
        )

        if not upstream_task_name:
            if self.consumes.get(requirement, True):
                raise ValueError(
                    f"Task '{self.task_name}' requires '{requirement.type.value}', "
                    f"but no upstream task was provided in 'depends_on'."
                )
            return ""  # Optional dependency not provided

        upstream_instance = self.parent.get_workflow_instance(upstream_task_name)
        if not any(p.satisfies(requirement) for p in upstream_instance.produces):
            raise TypeError("Upstream task does not satisfy required DataSignature")

        return f"output/{upstream_task_name}"

    def register_module(self, logical_name: str, default_module_name: str = ""):
        """
        Registers a dependency mapping for Snakemake 9.
        Paths are resolved using the explicit depends_on configuration.

        Arguments:
            logical_name: The logical name of the step depended on, e.g. "ngs_mapping" or "variant_calling"
            default_module_name: The default step to use if depends_on isn't specified explicitly,
                e.g., "somatic_variant_calling"
                (meaning that `depends_on: variant_calling: "somatic_variant_calling"` reads as "the variant_calling task is done by the task named somatic_variant_calling")
        """
        if not default_module_name:
            default_module_name = logical_name

        if self.depends_on is not None and hasattr(self.depends_on, logical_name):
            # Empty values in optional dependencies are treated as "not mapped".
            target_task_name = getattr(self.depends_on, logical_name) or default_module_name
        else:
            target_task_name = default_module_name

        if logical_name in self.modules:
            raise ValueError(f"Dependency {logical_name} already registered!")

        def resolve_dependency(path):
            # because replace_prefix only acts on intra-module files.
            if path.startswith("output/") or path.startswith("work/"):
                return f"{target_task_name}/" + path
            return os.path.join(target_task_name, path).replace("\\", "/")

        self.modules[logical_name] = resolve_dependency

    def get_upstream_path(self, logical_name: str, path: str, default_task_name: str = "") -> str:
        """Resolve *path* into the physical output directory of an upstream task.

        This is the inline alternative to the :meth:`register_module` +
        ``self.modules[logical_name](path)`` two-step pattern. Step parts can call this
        directly without first registering the module in ``__init__``:

        .. code-block:: python

            # old pattern (still works, no need to change immediately)
            self.register_module("ngs_mapping")
            ngs_mapping = self.parent.modules["ngs_mapping"]
            bam = ngs_mapping(f"output/{lib}/out/{lib}.bam")

            # new pattern (preferred for new code)
            bam = self.get_upstream_path("ngs_mapping", f"output/{lib}/out/{lib}.bam")

        Name resolution follows the same rules as :meth:`register_module`:

        1. If ``depends_on`` has an attribute named *logical_name*, use its value as the
           upstream task name.
        2. Otherwise fall back to *default_task_name* (or *logical_name* itself when
           *default_task_name* is empty).

        The task name is prepended to *path*, matching the Snakemake module namespace
        convention used everywhere in this codebase.

        Arguments:
            logical_name: The logical dependency name, e.g. ``"ngs_mapping"``.
            path: The local path within the upstream task, e.g.
                ``"output/lib/out/lib.bam"``.
            default_task_name: Fallback task name when *logical_name* is not found in
                ``depends_on``.  Defaults to *logical_name* itself.

        Returns:
            The namespaced path, e.g. ``"ngs_mapping/output/lib/out/lib.bam"``.
        """
        if not default_task_name:
            default_task_name = logical_name

        if self.depends_on is not None and hasattr(self.depends_on, logical_name):
            target_task_name = getattr(self.depends_on, logical_name) or default_task_name
        else:
            target_task_name = default_task_name

        if path.startswith("output/") or path.startswith("work/"):
            return f"{target_task_name}/{path}"
        return os.path.join(target_task_name, path).replace("\\", "/")

    def get_upstream_paths(
        self, req_field_name: str, **kwargs
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
                "ngs_mapping", library_name=wildcards.library_name
            )
            bam = alignments.bam
            bai = alignments.bai

        Arguments:
            req_field_name: The field name on ``self.config.depends_on`` that holds the upstream
                task reference. Must match a key in the ``depends_on`` Pydantic model.
            **kwargs: Forwarded verbatim to the upstream workflow's
                :py:meth:`get_output_paths` classmethod (e.g. ``library_name``, ``sample_name``).

        Returns:
            A ``pydantic.BaseModel`` instance when an expected-schema class is present in the
            ``Annotated`` metadata, otherwise a plain ``dict[str, str]``.

        Raises:
            AttributeError: If ``self.config`` has no ``depends_on`` attribute.
            ValueError: If the upstream task or its workflow class cannot be resolved.
        """
        if self.depends_on is None:
            raise AttributeError(
                f"Task '{self.task_name}' has no 'depends_on' configuration; "
                f"cannot resolve upstream field '{req_field_name}'."
            )

        target_task_name = getattr(self.depends_on, req_field_name)
        if not target_task_name:
            raise ValueError(
                f"Task '{self.task_name}': depends_on.{req_field_name} is empty or unset."
            )

        # Find the task entry to determine the step type
        target_task = next((t for t in self.w_config.tasks if t.name == target_task_name), None)
        step_type = target_task.step if target_task else target_task_name

        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        upstream_cls = WORKFLOW_REGISTRY.get(step_type)
        if upstream_cls is None:
            raise ValueError(
                f"Workflow class for step '{step_type}' (resolved from depends_on.{req_field_name}"
                f" = '{target_task_name}') not found in WORKFLOW_REGISTRY."
            )

        # Extract DataSignature and optional expected-schema from Annotated metadata
        field_info = self.config.depends_on.model_fields.get(req_field_name)
        required_sig: DataSignature | None = None
        expected_schema: type[pydantic.BaseModel] | None = None
        if field_info is not None:
            for meta in field_info.metadata:
                if isinstance(meta, DataSignature):
                    required_sig = meta
                elif isinstance(meta, ExpectedPathSchema):
                    expected_schema = meta.schema
                elif isinstance(meta, type) and issubclass(meta, pydantic.BaseModel):
                    expected_schema = meta

        # Delegate path construction to the upstream workflow classmethod
        local_paths = upstream_cls.get_output_paths(signature=required_sig, **kwargs)

        # Prepend upstream task name for Snakemake global namespace
        global_paths = {k: f"{target_task_name}/{v}" for k, v in local_paths.items()}

        if expected_schema is not None:
            return expected_schema(**global_paths)
        return global_paths

    def get_args(self, sub_step, action):
        """Return arguments for action of substep with given wildcards

        Delegates to the sub step object's get_args function
        """
        return self._get_sub_step(sub_step).get_args(action)

    def get_input_files(self, sub_step: str, action: str) -> Inputs | Callable[[Wildcards], Inputs]:
        """Return input files for action of substep with given wildcards

        Delegates to the sub step object's get_input_files function
        """
        input_files = self._get_sub_step(sub_step).get_input_files(action)
        if callable(input_files):

            def input_wrapper(*args, **kwargs):
                try:
                    return input_files(*args, **kwargs)
                except TypeError as e:
                    if kwargs and "unexpected keyword argument" in str(e):
                        if args:
                            return input_files(*args)
                        if "wildcards" in kwargs:
                            return input_files(kwargs["wildcards"])
                    raise

            return input_wrapper
        return input_files

    def get_output_files(self, sub_step: str, action: str) -> Outputs:
        """Return list of strings with output files/patterns

        Delegates to the sub step object's get_output_files function
        """
        return self._get_sub_step(sub_step).get_output_files(action)

    def get_params(self, sub_step: str, action: str) -> Any:
        """Return parameters

        Delegates to the sub step object's get_params function
        """
        return self.substep_dispatch(sub_step, "get_params", action)

    def get_resource(self, sub_step: str, action: str, resource_name: str) -> Any:
        """Get resource

        Delegates to the sub step object's get_resource function
        """
        return self.substep_dispatch(sub_step, "get_resource", action, resource_name)

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
        """Return path to the log file

        Delegates to the sub step object's get_log_file function
        """
        return self.substep_dispatch(sub_step, "get_log_file", action)

    def get_shell_cmd(self, sub_step: str, action: str, wildcards: Wildcards) -> str:
        """Return shell command for the pipeline sub step

        Delegates to the sub step object's get_shell_cmd function
        """
        return self.substep_dispatch(sub_step, "get_shell_cmd", action, wildcards)

    def run_locally(self, sub_step: str, action: str, wildcards: Wildcards) -> str:
        """Runs a function locally for the pipeline sub step

        Delegates to the sub step object's run_locally function
        """
        return self.substep_dispatch(sub_step, "run_locally", action, wildcards)

    def run(self, sub_step: str, action: str, wildcards: Wildcards) -> str:
        """Run command for the given action of the given sub step with the given wildcards

        Delegates to the sub step object's run function
        """
        return self._get_sub_step(sub_step).run(action, wildcards)

    def get_result_files(self) -> OutputFiles:
        """Return actual list of file names to build"""
        raise NotImplementedError("Implement me!")  # pragma: no cover

    def substep_getattr(self, step: str, name: str) -> Any:
        """Return attribute from substep"""
        return getattr(self._get_sub_step(step), name)

    def substep_dispatch(self, step: str, function: str, *args, **kwargs):
        """Dispatch call to function of sub step implementation"""
        return self.substep_getattr(step, function)(*args, **kwargs)

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
                self.config_lookup_paths,
                data_set.search_paths,
                data_set.search_patterns,
                data_set.type,
                data_set.is_background,
                data_set.naming_scheme,
                data_set.mixed_se_pe,
                data_set.sodar_uuid,
                data_set.sodar_title,
                data_set.pedigree_field,
            )

    def _load_data_search_infos(self) -> typing.Generator[DataSearchInfo, None, None]:
        """Use workflow and step configuration to yield ``DataSearchInfo`` objects"""
        for _, data_set in self.w_config.data_sets.items():
            yield DataSearchInfo(
                sheet_path=data_set.file,
                base_paths=self.config_lookup_paths,
                search_paths=self.config.search_paths,
                search_patterns=self.config.search_patterns,
                mixed_se_pe=False,
            )

    @classmethod
    def wrapper_path(cls, path: str) -> str:
        """Generate path to wrapper"""
        return "file://" + os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "snappy_wrappers",
                "wrappers",
                path,
            )
        )


class LinkInPathGenerator:
    """Helper class for generating paths to link in"""

    def __init__(
        self,
        work_dir,
        data_set_infos,
        config_paths,
        cache_file_name=".snappy_path_cache",
        preprocessed_path="",
    ):
        #: Working directory
        self.work_dir = work_dir
        #: Data set info list from configuration
        if preprocessed_path:
            self.data_set_infos = [
                self._update_datasetinfo(x, preprocessed_path) for x in data_set_infos
            ]
        else:
            self.data_set_infos = data_set_infos
        #: Path to configuration files, used for invalidating cache
        self.config_paths = config_paths
        #: Name of cache file to create
        self.cache_file_name = cache_file_name
        #: File system crawler to use
        invalidate_paths_list = self._merge_cache_invalidate_paths(self.data_set_infos)
        invalidate_paths_list += config_paths
        self.crawler = FileSystemCrawler(
            os.path.join(self.work_dir, self.cache_file_name), invalidate_paths_list
        )

    def run(
        self,
        folder_name,
        pattern_set_keys=("left", "right", "left_md5", "right_md5", "bam"),
    ):
        """Yield (src_path, path_infix, filename) one-by-one

        Cache is saved after the last iteration
        """
        # Iterate over data set infos and crawl file system
        filenames = set([])
        # TODO: crawling the actual data sheet of the current data set is enough!
        seen_root_paths = set()
        for info in self.data_set_infos:
            patterns = []
            # Build PatternSet objects, based on types in configuration
            for pat in info.search_patterns:
                if not isinstance(pat, (dict, MutableMapping)):
                    raise ValueError("search_patterns must be a dict!")  # pragma: no cover
                # Add patterns as found in configuration file: DataSetInfo.search_patterns
                patterns.append(PatternSet(pat.values(), names=pat.keys()))
                # Add MD5 files to search file as default
                pat_md5 = [pattern + ".md5" for pattern in pat.values()]
                pat_names_md5 = [pattern + "_md5" for pattern in pat.keys()]
                patterns.append(PatternSet(pat_md5, names=pat_names_md5))
            # Crawl all root paths, link in the resulting files
            for root_path in self._get_shell_cmd_root_paths(info):
                if root_path in seen_root_paths:
                    continue  # skip this root path
                seen_root_paths.add(root_path)
                for result in self.crawler.run(root_path, folder_name, patterns, info.mixed_se_pe):
                    res_dict = result.to_dict()
                    for key in pattern_set_keys:
                        if key not in res_dict:
                            continue  # skip if not found
                        path_infix = os.path.relpath(
                            os.path.dirname(res_dict[key]), result.base_folder
                        )
                        filename = os.path.basename(res_dict[key])
                        if res_dict[key] in filenames:
                            raise ValueError("Detected double link-in {}".format(filename))
                        filenames.add(filename)
                        src_dir = os.path.dirname(res_dict[key])
                        yield src_dir, path_infix, filename
        # Finally, save the cache
        self.crawler.save_cache()

    def _update_datasetinfo(self, data_set_info, preprocessed_path=""):
        return DataSetInfo(
            name=data_set_info.name,
            sheet_path=data_set_info.sheet_path,
            base_paths=data_set_info.base_paths,
            search_paths=[preprocessed_path] if preprocessed_path else data_set_info.search_paths,
            search_patterns=data_set_info.search_patterns,
            sheet_type=data_set_info.sheet_type,
            is_background=data_set_info.is_background,
            naming_scheme=data_set_info.naming_scheme,
            mixed_se_pe=data_set_info.mixed_se_pe,
            sodar_uuid=data_set_info.sodar_uuid,
            sodar_title=data_set_info.sodar_title,
        )

    @classmethod
    def _get_shell_cmd_root_paths(cls, info):
        for base_path in info.base_paths:
            if not os.path.exists(os.path.join(base_path, info.sheet_path)):
                continue  # skip this one
            for search_path in info.search_paths:
                yield os.path.abspath(
                    os.path.join(
                        os.path.dirname(os.path.join(base_path, info.sheet_path)),
                        search_path,
                    )
                )

    @classmethod
    def _merge_cache_invalidate_paths(cls, data_set_infos):
        """
        :param data_set_infos: List of DataSetInfo objects.
        :type data_set_infos: list

        :return: Returns list with paths that should be used to potentially
        invalidate a cache file based on the DataSetInfo. Method merges paths
        to a project tsv file as well as the search paths into a single list of strings.
        """
        # Initialise variable
        out_list = []
        # Iterate over DataSetInfo objects
        for info in data_set_infos:
            # Search paths - expects a list already
            out_list.extend(info.search_paths)

            # Sheet path
            # Only name of file is stored in config file (relative path used),
            # hence we need to find it in the base paths
            sheet_file_name = info.sheet_path  # expects a string
            base_paths = info.base_paths  # expects a list
            sheet_path = cls._find_sheet_file(sheet_file_name, base_paths)
            # Append if not None
            if sheet_path:
                out_list.append(sheet_path)

        # Return
        return out_list

    @classmethod
    def _find_sheet_file(cls, sheet_file_name, base_paths):
        """Method searches for sheet file in base paths.

        :param sheet_file_name: Sheet file name.
        :type sheet_file_name: str

        :param base_paths: List of strings with base paths.
        :type base_paths: list

        :return: Returns path to sheet file.
        """
        # Check if full path already
        if os.path.exists(sheet_file_name):
            return sheet_file_name
        # Iterate over base paths
        # Assumption: sheet file stored in the same level as config file,
        # i.e., one of the base paths.
        for base_p in base_paths:
            dir_path = os.path.realpath(base_p)
            # Find all files
            for item in os.listdir(dir_path):
                if sheet_file_name == item:
                    return os.path.join(dir_path, sheet_file_name)
        # If not found: None
        return None


def get_ngs_library_folder_name(sheets, library_name):
    """Return library's folder name

    The library is searched for based on the ``library_name``.  In the case of
    multiple NGS library matches, the first one is returned.
    """
    for sheet in sheets:
        try:
            ngs_library = sheet.crawl(sheet.name_generator.inverse(library_name))
        except SecondaryIDNotFoundException:
            continue  # skip, not in this sheet
        if ngs_library:
            try:
                return ngs_library.extra_infos["folderName"]
            except AttributeError:
                raise ValueError("No folderName extraInfos entry for {}".format(ngs_library.name))
    raise ValueError("Found no folders for NGS library of name {}".format(library_name))


class LinkInStepPart(BaseStepPart):
    """Link in the raw files, e.g. FASTQ files

    Depending on the configuration, the files are linked out after postprocessing
    """

    name = "link_in"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_pattern_out = "work/input_links/{library_name}/.done"
        self.preprocessed_path = self.parent.get_preprocessed_path()

        # Path generator.
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            cache_file_name=".snappy_path_cache",
            preprocessed_path=self.preprocessed_path,
        )

    def get_input_files(self, action):
        """Return required input files"""
        return []  # no input

    def get_output_files(self, action):
        assert action == "run", "Unsupported action"
        return touch(self.base_pattern_out)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking in the files"""
        assert action == "run", "Unsupported action"
        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""
        # Get base out path with the task prefix prepended
        out_path = os.path.dirname(task_prefix + self.base_pattern_out.format(**wildcards))
        # Get folder name of first library candidate
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.preprocessed_path:
            folder_name = wildcards.library_name
        # Perform the command generation
        lines = []
        tpl = (
            "mkdir -p {out_path}/{path_infix} && "
            "{{{{ test -h {out_path}/{path_infix}/{filename} || "
            "ln -sr {src_path}/{filename} {out_path}/{path_infix}; }}}}"
        )
        filenames = {}  # generated so far
        for src_path, path_infix, filename in self.path_gen.run(folder_name):
            new_path = os.path.join(out_path, path_infix, filename)
            if new_path in filenames:
                if filenames[new_path] == src_path:
                    continue  # ignore TODO: better correct this
                msg = "WARNING: Detected double output path {}"
                print(msg.format(filename), file=sys.stderr)
            filenames[new_path] = src_path
            lines.append(
                tpl.format(
                    src_path=src_path,
                    out_path=out_path,
                    path_infix=path_infix,
                    filename=filename,
                )
            )
        if not lines:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)
        return "\n".join(lines)

    def run_locally(self, action, wildcards):
        """Links fastq files"""
        assert action == "run", "Unsupported action"

        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""
        out_path = os.path.dirname(task_prefix + self.base_pattern_out.format(**wildcards))

        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.preprocessed_path:
            folder_name = wildcards.library_name

        filenames = self._create_all_symlinks(self.path_gen, folder_name, out_path)
        if not filenames:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)

    @staticmethod
    def _create_all_symlinks(
        path_generator,
        folder_name,
        out_path,
        pattern_set_keys=("left", "right", "left_md5", "right_md5", "bam"),
    ):
        filenames = {}  # generated so far
        for src_path, path_infix, filename in path_generator.run(folder_name, pattern_set_keys):
            new_path = os.path.join(out_path, path_infix, filename)
            if new_path in filenames:
                if filenames[new_path] == src_path:
                    continue  # ignore TODO: better correct this
                msg = "WARNING: Detected double output path {}"
                print(msg.format(filename), file=sys.stderr)
            filenames[new_path] = src_path
            # Create the symlink
            d = os.path.realpath(os.path.join(out_path, path_infix))
            link = os.path.join(d, filename)
            os.makedirs(d, exist_ok=True)
            if not os.path.islink(link):
                target = os.path.relpath(
                    os.path.join(os.path.realpath(src_path), filename), start=d
                )
                os.symlink(target, link)
        return filenames

    def run(self, action, wildcards):
        raise ImplementationUnavailableError(
            "run() not implemented for linking in reads"
        )  # pragma: no cover


class LinkInVcfExternalStepPart(LinkInStepPart):
    """Link in the external VCF files."""

    #: Step name
    name = "link_in_vcf_external"

    #: Class available actions
    actions = ("run",)

    #: Patterns set keys
    pattern_set_keys = ("vcf", "vcf_md5")

    def __init__(self, parent):
        super().__init__(parent)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking in the files

        The files are linked, keeping their relative paths to the item matching the "folderName"
        intact.
        """
        self._validate_action(action)
        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""
        # Define path generator
        path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_search_infos,
            self.parent.config_lookup_paths,
        )
        # Get base out path with the task prefix prepended
        out_path = os.path.dirname(task_prefix + self.base_pattern_out.format(**wildcards))
        # Perform the command generation
        lines = []
        tpl = (
            "mkdir -p {out_path}/{path_infix} && "
            "{{{{ test -h {out_path}/{path_infix}/{filename} || "
            "ln -sr {src_path}/{filename} {out_path}/{path_infix}; }}}}"
        )
        filenames = {}
        for src_path, path_infix, filename in path_gen.run(
            folder_name=wildcards.library_name, pattern_set_keys=self.pattern_set_keys
        ):
            new_path = os.path.join(out_path, path_infix, filename)
            if new_path in filenames:
                if filenames[new_path] == src_path:
                    continue  # ignore TODO: better correct this
                msg = "WARNING: Detected double output path {}"
                print(msg.format(filename), file=sys.stderr)
            filenames[new_path] = src_path
            lines.append(
                tpl.format(
                    src_path=src_path,
                    out_path=out_path,
                    path_infix=path_infix,
                    filename=filename,
                )
            )
        if not lines:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)
        return "\n".join(lines)

    def run_locally(self, action, wildcards):
        """Links fastq files

        The files are linked, keeping their relative paths to the item matching the "folderName"
        intact.
        """
        self._validate_action(action)
        task_prefix = f"{self.parent.task_name}/" if getattr(self.parent, "task_name", "") else ""
        # Define path generator
        path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_search_infos,
            self.parent.config_lookup_paths,
        )
        # Get base out path with the task prefix prepended
        out_path = os.path.dirname(task_prefix + self.base_pattern_out.format(**wildcards))
        filenames = self._create_all_symlinks(
            path_generator=path_gen,
            folder_name=wildcards.library_name,
            out_path=out_path,
            pattern_set_keys=self.pattern_set_keys,
        )
        if not filenames:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)


class LinkInBamExternalStepPart(LinkInVcfExternalStepPart):
    """Link in the external BAM files."""

    #: Step name
    name = "link_in_bam_external"

    #: Class available actions
    actions = ("run",)

    #: Patterns set keys
    pattern_set_keys = ("bam", "bam_md5")

    def __init__(self, parent):
        super().__init__(parent)
        self.base_pattern_out = "work/input_links/{library_name}/.done_bam_external"


class LinkInBaiExternalStepPart(LinkInVcfExternalStepPart):
    """Link in the external BAI files."""

    #: Step name
    name = "link_in_bai_external"

    #: Class available actions
    actions = ("run",)

    #: Patterns set keys
    pattern_set_keys = ("bai", "bai_md5")

    def __init__(self, parent):
        super().__init__(parent)
        self.base_pattern_out = "work/input_links/{library_name}/.done_bai_external"


class InputFilesStepPartMixin:
    """Mixin with predefined "get_input_files" function."""

    #: Whether to include path to PED file or not
    include_ped_file = None

    #: Class with input VCF file name
    prev_class = None

    #: Extensions of files to create as main payload
    ext_values = None

    #: Names of the files to create for the extension
    ext_names = None

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            if self.include_ped_file:
                yield (
                    "ped",
                    os.path.realpath(
                        "work/write_pedigree.{index_library}/out/{index_library}.ped"
                    ).format(**wildcards),
                )
            name_pattern = self.prev_class.name_pattern.replace(r",[^\.]+", "")
            tpl_path_out = os.path.join("work", name_pattern, "out", name_pattern)
            for key, ext in zip(self.ext_names, self.ext_values):
                yield key, tpl_path_out.format(ext=ext, **wildcards) + ext

        assert action == "run"
        return input_function
