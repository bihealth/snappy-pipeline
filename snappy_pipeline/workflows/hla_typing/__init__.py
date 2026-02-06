# -*- coding: utf-8 -*-
"""Implementation of the ``hla_typing`` step

The hla_typing step allows for the HLA typing from NGS read data (WGS, targeted DNA sequencing,
or RNA-seq).

==========
Step Input
==========

Gene fusion calling starts at the raw RNA-seq reads.  Thus, the input is very similar to one of
:ref:`ngs_mapping step <step_ngs_mapping>`.

See :ref:`ngs_mapping_step_input` for more information.

===========
Step Output
===========

HLA typing will be performed for all NGS libraries in all sample sheets.  For each combination
of HLA typer and library, a directory ``{hla_typer}.{lib_name}-{lib_pk}/out`` will be created.
Therein, the following files will be created:

- ``{hla_typer}.{lib_name}-{lib_pk}.txt``
- ``{hla_typer}.{lib_name}-{lib_pk}.txt.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- optitype.P001-N1-DNA1-WES1-4
    |   `-- out
    |       |-- optitype.P001-N1-DNA1-WES1-4.txt
    |       `-- optitype.P001-N1-DNA1-WES1-4.txt.md5
    [...]

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_hla_typing.rst

==========================
Available HLA Typing Tools
==========================

The following HLA typing tools are currently available

- ``"optitype"``
- ``"arcashla"``

"""

import os
import re
from collections import OrderedDict
from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand, Wildcards, InputFiles

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)

from .model import HlaTyping as HlaTypingConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".txt", ".txt.md5", ".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("txt", "txt_md5", "json", "json_md5")

#: HLA typing tools
HLA_TYPERS = ("optitype", "arcashla")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = HlaTypingConfigModel.default_config_yaml_string()


class OptiTypeStepPart(BaseStepPart):
    """HLA Typing using OptiType"""

    #: Step name
    name = "optitype"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/optitype.{{library_name}}/out/optitype.{{library_name}}{ext}"
        self.extensions = EXT_VALUES
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config.path_link_in,
        )

    @staticmethod
    def get_output_prefix():
        return ""

    @classmethod
    @dictify
    def get_input_files(cls, action):
        """Return input files"""
        assert action == "run"
        yield "done", "work/input_links/{library_name}/.done"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext)
        # add additional optitype output files
        for name, ext in {"tsv": "_result.tsv", "cov_pdf": "_coverage_plot.pdf"}.items():
            yield name, self.base_path_out.format(ext=ext)

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        self._validate_action(action)

        prefix = "work/{name}.{{library_name}}/log/{name}.{{library_name}}".format(name=self.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            result = {
                "input": {
                    "reads_left": list(
                        sorted(self._collect_reads(wildcards, wildcards.library_name, ""))
                    )
                },
                "seq_type": self._get_seq_type(wildcards),
            }
            reads_right = list(
                sorted(self._collect_reads(wildcards, wildcards.library_name, "right-"))
            )
            if reads_right:
                result["input"]["reads_right"] = reads_right
            result["use_discordant"] = "true" if self.config.optitype.use_discordant else "false"
            result["num_mapping_threads"] = self.config.optitype.num_mapping_threads
            result["max_reads"] = self.config.optitype.max_reads
            result["yara_error_rate"] = self.config.optitype.yara_mapper.error_rate
            result["yara_strata_rate"] = self.config.optitype.yara_mapper.strata_rate
            result["yara_sensitivity"] = self.config.optitype.yara_mapper.sensitivity
            return result

        assert action == "run", "Unsupported actions"
        return args_function

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config.path_link_in:
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)

    def _get_seq_type(self, wildcards):
        """Return sequence type for the library name in wildcards"""
        library = self.parent.ngs_library_name_to_ngs_library[wildcards.library_name]
        return library.test_sample.extra_infos.get("extractionType", "DNA").lower()

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        return ResourceUsage(
            threads=6,
            time="40:00:00",  # 40 hours
            memory="45000M",
        )


class ArcasHlaStepPart(BaseStepPart):
    """HLA Typing using arcasHLA"""

    #: Step name
    name = "arcashla"

    #: Class available actions
    actions = ("run",)

    PAIRED_PATTERN: re.Pattern = re.compile(
        r"^(?P<passPairs>[0-9]+) \+ (?P<failPairs>[0-9]+) paired in sequencing$"
    )

    def __init__(self, parent):
        super().__init__(parent)
        self.mapper = self.config.arcashla.mapper
        self.base_path_out = (
            "work/{mapper}.{name}.{{library_name}}/out/{mapper}.{name}.{{library_name}}{ext}"
        )
        self.extensions = EXT_VALUES

    def get_input_files(self, action):
        """Return input files"""

        @dictify
        def input_function(wildcards):
            yield "ref_done", "work/arcashla.prepare_reference/out/.done"
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            tpl = f"{self.mapper}.{wildcards.library_name}"
            yield "bam", ngs_mapping(f"output/{tpl}/out/{tpl}.bam")
            yield "flagstats", ngs_mapping(f"output/{tpl}/report/bam_qc/{tpl}.bam.flagstats.txt")

        assert action == "run"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        self._validate_action(action)
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext, mapper=self.mapper, name=self.name)

    def get_output_prefix(self):
        return "%s." % self.mapper

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        args = dict(self.config.arcashla.model_dump(by_alias=True))
        del args["mapper"]
        if args["drop_iterations"] is None:
            del args["drop_iterations"]
        if args["avg"] is None:
            del args["avg"]
        if args["std"] is None:
            del args["std"]

        with open(input["flagstats"], "rt") as f:
            for line in f:
                m = self.PAIRED_PATTERN.match(line.strip())
                if m:
                    try:
                        n = int(m.group("passPairs")) + int(m.group("failPairs"))
                    except ValueError:
                        n = 0
                    args["single"] = n == 0
                    break

        return args

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        self._validate_action(action)

        prefix = (
            "work/{mapper}.{name}.{{library_name}}/log/{mapper}.{name}.{{library_name}}".format(
                mapper=self.mapper, name=self.name
            )
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        return ResourceUsage(
            threads=4,
            time="60:00:00",  # 60 hours
            memory="15000M",
        )


class HlaLaStepPart(BaseStepPart):
    """HLA Typing using HLA-LA"""

    #: Step name
    name = "hla_la"

    #: Class available actions
    actions = ("prepare_graph", "prepare_reference", "run")

    FASTA_PATTERN: re.Pattern = re.compile(r"\.fa(sta)?(\.gz)?$")

    def __init__(self, parent):
        super().__init__(parent)
        self.mapper = self.config.hla_la.mapper
        self.base_path_out = (
            "work/{mapper}.{name}.{{library_name}}/out/{mapper}.{name}.{{library_name}}{ext}"
        )
        self.extensions = EXT_VALUES

        if self.config.hla_la.path_graph:
            self.path_graph = self.config.hla_la.path_graph
        else:
            self.path_graph = "work/hla_la.prepareGraph/out/.done"
        self.path_reference = os.path.join(
            os.path.dirname(self.path_graph),
            "knownReferences",
            self.FASTA_PATTERN.sub(
                ".txt",
                os.path.basename(self.w_config.static_data_config.reference.path),
            ),
        )

    def get_input_files(self, action):
        """Return input files"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_prepare_reference(self, wildcards: Wildcards):
        yield "path_graph", self.path_graph
        yield "reference", self.w_config.static_data_config.reference.path + ".fai"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "path_graph", self.path_graph
        yield "reference", self.path_reference
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = f"{self.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{tpl}/out/{tpl}.bam")

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        match action:
            case "prepare_graph":
                yield "done", self.path_graph
            case "prepare_reference":
                yield "reference", self.path_reference
            case "run":
                for name, ext in zip(EXT_NAMES, EXT_VALUES):
                    yield (
                        name,
                        self.base_path_out.format(ext=ext, mapper=self.mapper, name=self.name),
                    )
            case _:
                raise UnsupportedActionException(
                    f"Unsupported action {action} for tool {self.name}"
                )

    def get_output_prefix(self):
        return "%s." % self.mapper

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_prepare_reference(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"start": self.config.hla_la.start, "end": self.config.hla_la.end}

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"sample_id": wildcards.library_name}

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        self._validate_action(action)

        prefix = (
            "work/{mapper}.{name}.{{library_name}}/log/{mapper}.{name}.{{library_name}}".format(
                mapper=self.mapper, name=self.name
            )
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="60:00:00",  # 60 hours
            memory="60000M",
        )


class HlaTypingWorkflow(BaseStep):
    """Perform HLA Typing"""

    #: Step name
    name = "hla_typing"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, config_model_class=HlaTypingConfigModel)
        sub_steps = [LinkInStep, LinkOutStepPart, OptiTypeStepPart, ArcasHlaStepPart, HlaLaStepPart]
        self.register_sub_step_classes(tuple(sub_steps))
        #: Mapping from library name to library object
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library
        # Register sub workflows
        for tool in self.config.tools.dna + self.config.tools.rna:
            if self.config.get(tool).get("mapper", None):
                self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
                break

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        from os.path import join

        name_pattern = "{prefix}{hla_typer}.{ngs_library.name}"
        yield from self._yield_result_files(
            join("output", name_pattern, "out", name_pattern + "{ext}"), ext=EXT_VALUES
        )
        log_ext = [e + m for e in ("log", "conda_list.txt", "conda_info.txt") for m in ("", ".md5")]
        yield from self._yield_result_files(
            join("output", name_pattern, "log", name_pattern + ".{ext}"), ext=log_ext
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                extraction_type = ngs_library.test_sample.extra_infos.get(
                    "extractionType", "DNA"
                ).lower()
                for tool in self.config.tools.get(extraction_type):
                    yield from expand(
                        tpl,
                        prefix=self.sub_steps[tool].get_output_prefix(),
                        hla_typer=[tool],
                        ngs_library=[ngs_library],
                        **kwargs,
                    )
