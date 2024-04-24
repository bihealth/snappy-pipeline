# -*- coding: utf-8 -*-
"""Base classes for the actual pipeline steps"""

import contextlib
import datetime
import itertools
import os
import os.path
import sys
import tempfile
import typing
from collections import OrderedDict
from collections.abc import MutableMapping
from fnmatch import fnmatch
from functools import lru_cache
from io import StringIO

import attr
import ruamel.yaml as ruamel_yaml
from biomedsheets import io_tsv
from biomedsheets.io import SheetBuilder, json_loads_ordered
from biomedsheets.models import SecondaryIDNotFoundException
from biomedsheets.naming import NAMING_SCHEMES, NAMING_SECONDARY_ID_PK, name_generator_for_scheme
from biomedsheets.ref_resolver import RefResolver
from biomedsheets.shortcuts import (
    donor_has_dna_ngs_library,
    write_pedigree_to_ped,
    write_pedigrees_to_ped,
)
from snakemake.io import touch

from snappy_pipeline.base import (
    MissingConfiguration,
    UnsupportedActionException,
    merge_dicts,
    merge_kwargs,
    print_config,
    print_sample_sheets,
)
from snappy_pipeline.find_file import FileSystemCrawler, PatternSet
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract.pedigree import append_pedigree_to_ped
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


class BaseStepPart:
    """Base class for a part of a pipeline step"""

    name = "<base step>"

    #: The actions available in the class.
    actions: typing.Tuple[str] = None

    #: Default resource usage for actions that are not given in ``resource_usage``.
    default_resource_usage: ResourceUsage = ResourceUsage(
        threads=1,
        time="01:00:00",
        memory="2G",  # 1h
    )

    #: Configure resource usage here that should not use the default resource usage from
    #: ``default_resource_usage``.
    resource_usage: typing.Dict[str, ResourceUsage] = {}

    def __init__(self, parent):
        self.name = self.__class__.name
        self.parent = parent
        self.config = parent.config
        self.w_config = parent.w_config

    def _validate_action(self, action):
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
        """Return the resource usage for the given action."""
        if action not in self.actions:
            raise ValueError(f"Invalid {action} not in {self.actions}")
        return self.resource_usage.get(action, self.default_resource_usage)

    @staticmethod
    def get_default_partition() -> str:
        """Helper that returns the default partition."""
        return os.getenv("SNAPPY_PIPELINE_PARTITION")

    def get_resource(self, action: str, resource_name: str):
        """Return the amount of resources to be allocated for the given action.

        :param action: The action to return the resource requirement for.
        :param resource_name: The name to return the resource for.
        """
        if resource_name not in ("threads", "time", "memory", "partition", "tmpdir"):
            raise ValueError(f"Invalid resource name: {resource_name}")
        resource_usage = self.get_resource_usage(action)
        if resource_name == "tmpdir" and not resource_usage.tmpdir:
            return self.parent.get_tmpdir()
        if resource_name == "partition" and not resource_usage.partition:
            return self.get_default_partition()
        else:
            return getattr(resource_usage, resource_name)

    def get_args(self, action):
        """Return args for the given action of the sub step"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_input_files(self, action):
        """Return input files for the given action of the sub step"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_output_files(self, action):
        """Return output files for the given action of the sub step and"""
        raise NotImplementedError("Called abstract method. Override me!")  # pragma: no cover

    def get_log_file(self, action):
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

    def get_shell_cmd(self, action, wildcards):  # NOSONAR
        """Return shell command for the given action of the sub step and the given wildcards"""
        raise ImplementationUnavailableError(
            "Override this method before calling it!"
        )  # pragma: no cover

    def run(self, action, wildcards):  # NOSONAR
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

    def __init__(self, parent, require_dna_ngs_library=False, only_trios=False):
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
                                in_trio |= {donor.name, donor.father.name, donor.mother.name}
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
            # if "ngs_mapping" not in self.parent.sub_workflows:
            #     return  # early exit
            # Get shortcut to NGS mapping sub workflow
            # Get names of primary libraries of the selected pedigree.  The pedigree is selected
            # by the primary DNA NGS library of the index.
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            if not pedigree.index or not pedigree.index.dna_ngs_library:
                msg = "INFO: pedigree without index (names: {})"  # pragma: no cover
                donor_names = list(sorted(d.name for d in pedigree.donors))
                print(msg.format(donor_names), file=sys.stderr)  # pragma: no cover
                return
            mappers = self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"]
            tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
            for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
                library_name = donor.dna_ngs_library.name
                for mapper in mappers:
                    path = tpl.format(
                        library_name=library_name, mapper=mapper, ext=".bam", **wildcards
                    )
                    yield Path("../ngs_mapping") / path

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

    def run(self, wildcards, output):
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
                self.index_ngs_library_to_pedigree[wildcards.index_ngs_library], str(output)
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

        def input_function(wildcards):
            """Helper wrapper function"""
            result = self.base_path_in.format(**wildcards)
            for pattern in self.disable_patterns:
                if fnmatch(result, pattern):
                    raise ValueError("Blocking match...")
            return result

        assert action == "run", "Unsupported action"
        return input_function

    def get_output_files(self, action):
        """Return output file pattern"""
        assert action == "run", "Unsupported action"
        return self.base_pattern_out

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out"""
        assert action == "run", "Unsupported action"
        tpl = "test -h {out} || ln -sr {in_} {out}"
        in_ = self.base_path_in.replace("{", "{wildcards.")
        out = self.base_path_out.replace("{", "{wildcards.")
        return tpl.format(in_=in_, out=out)


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
        sheet.json_data["extraInfoDefs"]["is_background"] = {"type": "boolean", "default": False}
        sheet.extra_infos["is_background"] = flag
        return sheet


@attr.s(frozen=True, auto_attribs=True)
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
    name = None

    #: Override with the sheet shortcut class to use
    sheet_shortcut_class = None

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

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        work_dir,
        previous_steps=None,
    ):
        self.name = self.__class__.name
        #: Tuple with absolute paths to configuration files read
        self.config_paths = config_paths
        #: Absolute path to directory of where to perform work
        self.work_dir = work_dir
        #: Classes of previously executed steps, used for merging their default configuration as
        #: well.
        self.previous_steps = tuple(previous_steps or [])
        #: Snakefile "workflow" object
        self.workflow = workflow
        self.modules = {}
        #: Merge default configuration with true configuration
        self.w_config = config
        self.w_config.update(self._update_config(config))
        self.config = self.w_config["step_config"].get(self.name, OrderedDict())
        #: Paths with configuration paths, important for later retrieving sample sheet files
        self.config_lookup_paths = config_lookup_paths
        self.sub_steps = {}
        self.data_set_infos = list(self._load_data_set_infos())
        # Check configuration
        self._check_config()
        #: Shortcut to the BioMed SampleSheet objects
        self.sheets = [info.sheet for info in self.data_set_infos]
        #: Shortcut BioMed SampleSheet keyword arguments
        sheet_kwargs_list = [
            merge_kwargs(
                first_kwargs=self.sheet_shortcut_kwargs, second_kwargs=info.pedigree_field_kwargs
            )
            for info in self.data_set_infos
        ]
        #: Shortcut sheets
        self.shortcut_sheets = []
        klass = self.__class__.sheet_shortcut_class
        for sheet, kwargs in zip(self.sheets, sheet_kwargs_list):
            kwargs = kwargs or {}
            kwargs = {k: v for k, v in kwargs.items() if k in klass.supported_kwargs}
            self.shortcut_sheets.append(
                klass(sheet, *(self.__class__.sheet_shortcut_args or []), **kwargs)
            )
        # Setup onstart/onerror/onsuccess hooks
        self._setup_hooks()

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

    def _update_config(self, config):
        """Update configuration config with the configuration returned by subclass'
        ``default_config_yaml()`` and return
        """
        result = OrderedDict()
        for cls in itertools.chain([self.__class__], self.previous_steps):
            result = merge_dicts(
                result, _cached_yaml_round_trip_load_str(cls.default_config_yaml())
            )
        return merge_dicts(result, config)

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
            if entry in handle:
                handle = handle[entry]
                so_far.append(entry)
            else:
                tpl = 'Missing configuration ("{full_path}", got up to "{so_far}"): {msg}'.format(
                    full_path="/".join(config_keys), so_far="/".join(so_far), msg=msg
                )
                raise e_class(tpl)

    def register_sub_step_classes(self, classes):
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
            obj.check_config()
            self.sub_steps[klass.name] = obj

    def register_module(self, step_name: str, prefix: os.PathLike, module_name: str | None = None):
        """
        Register workflow with given pipeline ``step_name``, using the given ``prefix``.
        This requires importing the respective workflow in the Snakefile
        (since the module API is not intended to be used programmatically).
        For example:

        ```
        module ngs_mapping:
            snakefile:
                "../ngs_mapping/Snakefile"
            config:
                wf.w_config
            prefix:
                wf.w_config["step_config"]["your_workflow"].get("path_ngs_mapping", "../ngs_mapping")


        use rule * from ngs_mapping
        ```

        Optionally, the module name can be given separate from ``step_name`` (the default)
        value for it.
        """
        module_name = module_name or step_name
        if module_name in self.modules:
            raise ValueError("Sub workflow {} already registered!".format(module_name))
        self.modules[module_name] = lambda path: os.path.join(prefix, path)

    def get_args(self, sub_step, action):
        """Return arguments for action of substep with given wildcards

        Delegates to the sub step object's get_input_files function
        """
        return self._get_sub_step(sub_step).get_args(action)

    def get_input_files(self, sub_step, action):
        """Return input files for action of substep with given wildcards

        Delegates to the sub step object's get_input_files function
        """
        return self._get_sub_step(sub_step).get_input_files(action)

    def get_output_files(self, sub_step, action):
        """Return list of strings with output files/patterns

        Delegates to the sub step object's get_output_files function
        """
        return self._get_sub_step(sub_step).get_output_files(action)

    def get_params(self, sub_step, action):
        """Return parameters

        Delegates to the sub step object's get_params function
        """
        return self.substep_dispatch(sub_step, "get_params", action)

    def get_resource(self, sub_step, action, resource_name):
        """Get resource

        Delegates to the sub step object's get_resource function
        """
        return self.substep_dispatch(sub_step, "get_resource", action, resource_name)

    def get_tmpdir(self):
        """Return temporary directory.

        To be used directly or via get_resource("step", "action", "tmpdir")

        1. Try to evaluate global_config/tmpdir. Interpret $-variables from environment.
           Provides the current date as $TODAY.
        2. If this fails, try to use environment variable TMPDIR.
        3. If this fails, use tempfile.gettempdir(), same as Snakemake default.
        """
        tmpdir = self.w_config.get("global_config", {}).get("tmpdir", None)
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

    def get_log_file(self, sub_step, action):
        """Return path to the log file

        Delegates to the sub step object's get_log_file function
        """
        return self.substep_dispatch(sub_step, "get_log_file", action)

    def get_shell_cmd(self, sub_step, action, wildcards):
        """Return shell command for the pipeline sub step

        Delegates to the sub step object's get_shell_cmd function
        """
        return self.substep_dispatch(sub_step, "get_shell_cmd", action, wildcards)

    def run(self, sub_step, action, wildcards):
        """Run command for the given action of the given sub step with the given wildcards

        Delegates to the sub step object's run function
        """
        return self._get_sub_step(sub_step).get_shell_cmd(action, wildcards)

    def get_result_files(self):
        """Return actual list of file names to build"""
        raise NotImplementedError("Implement me!")  # pragma: no cover

    def substep_getattr(self, step, name):
        """Return attribute from substep"""
        return getattr(self._get_sub_step(step), name)

    def substep_dispatch(self, step, function, *args, **kwargs):
        """Dispatch call to function of sub step implementation"""
        return self.substep_getattr(step, function)(*args, **kwargs)

    def _get_sub_step(self, sub_step):
        if sub_step in self.sub_steps:
            return self.sub_steps[sub_step]
        else:
            raise ValueError(
                'Could not find sub step "{}" in workflow step "{}"'.format(sub_step, self.name)
            )  # pragma: no cover

    def _load_data_set_infos(self):
        """Load BioMed Sample Sheets as given by configuration and yield them"""
        for name, data_set in self.w_config["data_sets"].items():
            yield DataSetInfo(
                name,
                data_set["file"],
                self.config_lookup_paths,
                data_set["search_paths"],
                data_set["search_patterns"],
                data_set["type"],
                data_set.get("is_background", False),
                data_set.get("naming_scheme", NAMING_SECONDARY_ID_PK),
                data_set.get("mixed_se_pe", False),
                data_set.get("sodar_uuid", None),
                data_set.get("sodar_title", None),
                data_set.get("pedigree_field", None),
            )

    def _load_data_search_infos(self):
        """Use workflow and step configuration to yield ``DataSearchInfo`` objects"""
        for _, data_set in self.w_config["data_sets"].items():
            yield DataSearchInfo(
                sheet_path=data_set["file"],
                base_paths=self.config_lookup_paths,
                search_paths=self.config["search_paths"],
                search_patterns=self.config["search_patterns"],
                mixed_se_pe=False,
            )

    @classmethod
    def wrapper_path(cls, path):
        """Generate path to wrapper"""
        return "file://" + os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), "..", "..", "..", "snappy_wrappers", "wrappers", path
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

    def run(self, folder_name, pattern_set_keys=("left", "right", "left_md5", "right_md5", "bam")):
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
                        os.path.dirname(os.path.join(base_path, info.sheet_path)), search_path
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


# TODO: Rename to LinkInStepPart
class LinkInStep(BaseStepPart):
    """Link in the raw files, e.g. FASTQ files

    Depending on the configuration, the files are linked out after postprocessing
    """

    name = "link_in"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_pattern_out = "work/input_links/{library_name}/.done"

        # The key 'path_link_in' is only defined for pipelines that could used preprocessed
        # FASTQ files. That doesn't make sense for pipelines that are using externally generated
        # data already.
        try:
            preprocessed_path = self.config["path_link_in"]
        except KeyError:
            preprocessed_path = ""

        # Path generator.
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            cache_file_name=".snappy_path_cache",
            preprocessed_path=preprocessed_path,
        )

    def get_input_files(self, action):
        """Return required input files"""
        return []  # no input

    def get_output_files(self, action):
        assert action == "run", "Unsupported action"
        return touch(self.base_pattern_out)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking in the files

        The files are linked, keeping their relative paths to the item matching the "folderName"
        intact.
        """
        assert action == "run", "Unsupported action"
        # Get base out path
        out_path = os.path.dirname(self.base_pattern_out.format(**wildcards))
        # Get folder name of first library candidate
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config["path_link_in"]:
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
                    src_path=src_path, out_path=out_path, path_infix=path_infix, filename=filename
                )
            )
        if not lines:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)
        return "\n".join(lines)

    def run(self, action, wildcards):
        raise ImplementationUnavailableError(
            "run() not implemented for linking in reads"
        )  # pragma: no cover


class LinkInVcfExternalStepPart(LinkInStep):
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
        # Define path generator
        path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_search_infos, self.parent.config_lookup_paths
        )
        # Get base out path
        out_path = os.path.dirname(self.base_pattern_out.format(**wildcards))
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
                    src_path=src_path, out_path=out_path, path_infix=path_infix, filename=filename
                )
            )
        if not lines:
            msg = "Found no files to link in for {}".format(dict(**wildcards))
            print(msg, file=sys.stderr)
            raise Exception(msg)
        return "\n".join(lines)


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
