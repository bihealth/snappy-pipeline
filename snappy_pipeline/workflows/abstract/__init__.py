# -*- coding: utf-8 -*-
"""Base classes for the actual pipeline steps"""

from collections import OrderedDict
from collections.abc import MutableMapping
from fnmatch import fnmatch
from functools import lru_cache
from io import StringIO
import itertools
import os
import os.path
import sys

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
import ruamel.yaml as yaml
from snakemake.io import touch

from snappy_pipeline.base import (
    MissingConfiguration,
    merge_dicts,
    merge_kwargs,
    print_config,
    print_sample_sheets,
    snakefile_path,
)
from snappy_pipeline.find_file import FileSystemCrawler, PatternSet
from snappy_pipeline.utils import dictify, listify

#: String constant with bash command for redirecting stderr to ``{log}`` file
STDERR_TO_LOG_FILE = r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file and enable printing executed commands
exec 2> >(tee -a "{log}")
set -x
# -----------------------------------------------------------------------------

""".lstrip()


class ImplementationUnavailableError(NotImplementedError):
    """Raised when a function that is to be overridden optionally is called

    This is provided as an alternative to ``NotImplementedError`` as the Python linters warn if
    a class does not override functions throwing ``NotImplementedError``.
    """


class BaseStepPart:
    """Base class for a part of a pipeline step"""

    name = "<base step>"

    def __init__(self, parent):
        self.name = self.__class__.name
        self.parent = parent
        self.config = parent.config
        self.w_config = parent.w_config

    def update_cluster_config(self, cluster_config):
        """Override and configure the cluster resource requirements for all rules that this
        pipeline step part uses
        """

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

    name = "write_pedigree"

    def __init__(self, parent, require_dna_ngs_library=False, only_trios=False):
        super().__init__(parent)
        #: Whether or not to prevent writing out of samples with out NGS library.
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

        @listify
        def get_input_files(wildcards):
            if "ngs_mapping" not in self.parent.sub_workflows:
                return  # early exit
            # Get shorcut to NGS mapping sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
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
                    yield ngs_mapping(path)

        assert action == "run", "Unsupported actions"
        return get_input_files

    def get_output_files(self, action):
        assert action == "run"
        return "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        # return os.path.realpath(
        #    "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        # )

    def run(self, wildcards, output):
        """Write out the pedigree information"""
        if wildcards.index_ngs_library == "whole_cohort":
            write_pedigrees_to_ped(self.index_ngs_library_to_pedigree.values(), str(output))
        else:
            write_pedigree_to_ped(
                self.index_ngs_library_to_pedigree[wildcards.index_ngs_library], str(output)
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
    return yaml.round_trip_load(StringIO(str_value))


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

        :param base_paths:  All base paths of all configuration, too look for ``sheet_path``.
        :type base_paths: list

        :param search_paths: Search paths for the files in the sample sheet.

        :param search_patterns: Search patterns. Example: "{left: '**/*_R1_*.fastq.gz',
        right: '**/*_R2_*.fastq.gz'}".
        :type search_patterns: dict

        :param sheet_type: Explicite sheet type (e.g. "matched_cancer"), if any.  Otherwise, will
        attempt to load from sheet.
        :type sheet_type: str

        :param is_background: Whether or not the data set info is to be used only for background.
        :type is_background: bool

        :param naming_scheme: Selected naming schema: either 'secondary_id_pk' or
        'only_secondary_id'.
        :type naming_scheme: str

        :param mixed_se_pe: Whether or not mixing SE and PE data sets is allowed.

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
        #: All base paths of all configuration, too look for ``sheet_path``
        self.base_paths = base_paths
        #: Search paths for the files in the sample sheet
        self.search_paths = list(search_paths)
        #: Search patterns
        self.search_patterns = search_patterns
        #: Explicite sheet type (e.g. "matched_cancer"), if any.  Otherwise, will attempt to load
        # from sheet.
        self.sheet_type = sheet_type
        #: Whether or not the data set info is to be used only for background
        self.is_background = is_background
        #: Selected naming schema
        if naming_scheme not in NAMING_SCHEMES:
            raise ValueError("Invalid naming scheme: {}".format(naming_scheme))  # pragma: no cover
        self.naming_scheme = naming_scheme
        #: Whether or not mixing SE and PE data sets is allowed.
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
        cluster_config,
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
        #: Merge default configuration with true configuration
        self.w_config = config
        self.w_config.update(self._update_config(config))
        self.config = self.w_config["step_config"].get(self.name, OrderedDict())
        #: Cluster configuration dict
        self.cluster_config = cluster_config
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
        #: Functions from sub workflows, can be used to generate output paths into these workflows
        self.sub_workflows = {}

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
        :type config_keys: list

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

    def update_cluster_config(self):
        """Update cluster configuration for rule "__default__"

        The sub parts' ``update_cluster_config()`` routines are called on creation in
        ``register_sub_step_classes()``.
        """
        self.cluster_config["__default__"] = {"mem": 4 * 1024, "time": "12:00", "ntasks": 1}

    def register_sub_step_classes(self, classes):
        """Register an iterable of sub step classes

        Initializes objects in ``self.sub_steps`` dict
        """
        self.update_cluster_config()
        for pair_or_class in classes:
            try:
                klass, args = pair_or_class
            except TypeError:
                klass = pair_or_class
                args = ()
            obj = klass(self, *args)
            obj.update_cluster_config(self.cluster_config)
            obj.check_config()
            self.sub_steps[klass.name] = obj

    def register_sub_workflow(self, step_name, workdir, sub_workflow_name=None):
        """Register workflow with given pipeline ``step_name`` and in the given ``workdir``.

        Optionally, the sub workflow name can be given separate from ``step_name`` (the default)
        value for it.
        """
        sub_workflow_name = sub_workflow_name or step_name
        if sub_workflow_name in self.sub_workflows:
            raise ValueError("Sub workflow {} already registered!".format(sub_workflow_name))
        if os.path.isabs(workdir):
            abs_workdir = workdir
        else:
            abs_workdir = os.path.realpath(os.path.join(os.getcwd(), workdir))
        self.workflow.subworkflow(
            sub_workflow_name,
            workdir=abs_workdir,
            snakefile=snakefile_path(step_name),
            configfile=abs_workdir + "/" + "config.yaml",
        )
        self.sub_workflows[sub_workflow_name] = self.workflow.globals[sub_workflow_name]

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
        self, work_dir, data_set_infos, config_paths, cache_file_name=".snappy_path_cache"
    ):
        #: Working directory
        self.work_dir = work_dir
        #: Data set info list from configuration
        self.data_set_infos = data_set_infos
        #: Path to configuration files, used for invalidating cache
        self.config_paths = config_paths
        #: Name of cache file to create
        self.cache_file_name = cache_file_name
        #: File system crawler to use
        invalidate_paths_list = self._merge_cache_invalidate_paths(data_set_infos)
        invalidate_paths_list += config_paths
        self.crawler = FileSystemCrawler(
            os.path.join(self.work_dir, self.cache_file_name), invalidate_paths_list
        )

    def run(self, folder_name, pattern_set_keys=("left", "right", "bam")):
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
                patterns.append(PatternSet(pat.values(), names=pat.keys()))
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
            out_list.extend(getattr(info, "search_paths"))

            # Sheet path
            # Only name of file is stored in config file (relative path used),
            # hence we need to find it in the base paths
            sheet_file_name = getattr(info, "sheet_path")  # expects a string
            base_paths = getattr(info, "base_paths")  # expects a list
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
        # Path generator.
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_set_infos, self.parent.config_lookup_paths
        )

    def get_input_files(self, action):
        """Return required input files"""
        return []  # no input

    def get_output_files(self, action):
        """Return output files that are generated by snappy-gatk_post_bam"""
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
                yield "ped", os.path.realpath(
                    "work/write_pedigree.{index_library}/out/{index_library}.ped"
                ).format(**wildcards)
            name_pattern = self.prev_class.name_pattern.replace(r",[^\.]+", "")
            tpl_path_out = os.path.join("work", name_pattern, "out", name_pattern)
            for key, ext in zip(self.ext_names, self.ext_values):
                yield key, tpl_path_out.format(ext=ext, **wildcards) + ext

        assert action == "run"
        return input_function
