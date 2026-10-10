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

HLA typing will be performed for all NGS libraries in all sample sheets. For each
library, a directory ``{lib_name}-{lib_pk}/out`` will be created.
Therein, the following files will be created:

- ``{lib_name}-{lib_pk}.txt``
- ``{lib_name}-{lib_pk}.txt.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- P001-N1-DNA1-WES1-4
    |   `-- out
    |       |-- P001-N1-DNA1-WES1-4.txt
    |       `-- P001-N1-DNA1-WES1-4.txt.md5
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

import re
from collections import OrderedDict

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import HlaTyping as HlaTypingConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".txt", ".txt.md5", ".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("txt", "txt_md5", "json", "json_md5")

#: HLA typing tools
HLA_TYPERS = ("optitype", "arcashla")

#: Default configuration for the hla_typing schema


class OptiTypeStepPart(BaseStepPart):
    """HLA Typing using OptiType"""

    #: Step name
    name = "optitype"

    supported_extraction_types = ("dna", "rna")

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"
        self.extensions = EXT_VALUES

    @staticmethod
    def get_output_prefix():
        return ""

    def _get_input_files_run(self, wildcards):
        """Return input files"""
        return reads_input_files(self.parent, wildcards.library_name)

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        if self.name != str(self.config.tool):
            return {}
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext)
        # add additional optitype output files
        for name, ext in {"tsv": ".result.tsv", "cov_pdf": ".coverage_plot.pdf"}.items():
            yield name, self.base_path_out.format(ext=ext)

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        self._validate_action(action)

        if self.name != str(self.config.tool):
            return {}
        prefix = "work/{library_name}/log/{library_name}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result["seq_type"] = self._get_seq_type(wildcards)
        result["use_discordant"] = "true" if self.config.optitype.use_discordant else "false"
        result["num_mapping_threads"] = self.config.optitype.num_mapping_threads
        result["max_reads"] = self.config.optitype.max_reads
        result["yara_error_rate"] = self.config.optitype.yara_mapper.error_rate
        result["yara_strata_rate"] = self.config.optitype.yara_mapper.strata_rate
        result["yara_sensitivity"] = self.config.optitype.yara_mapper.sensitivity
        return result

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
            runtime="40h",  # 40 hours
            mem="45000MB",
        )


class ArcasHlaStepPart(BaseStepPart):
    """HLA Typing using arcasHLA"""

    #: Step name
    name = "arcashla"

    supported_extraction_types = ("rna",)

    #: Class available actions
    actions = ("run",)

    PAIRED_PATTERN: re.Pattern = re.compile(
        r"^(?P<passPairs>[0-9]+) \+ (?P<failPairs>[0-9]+) paired in sequencing$"
    )

    def __init__(self, parent):
        super().__init__(parent)
        self.mapper = self.config.arcashla.mapper
        self.base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"
        self.extensions = EXT_VALUES

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        yield "ref_done", "work/prepare_reference/out/.done"
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        if self.name != str(self.config.tool):
            return {}
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext, mapper=self.config.arcashla.mapper)

    def get_output_prefix(self):
        return ""

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        self._validate_action(action)
        prefix = "work/{library_name}/log/{library_name}"
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
            runtime="60h",  # 60 hours
            mem="15000MB",
        )


class HlaTypingWorkflow(BaseStep):
    """Perform HLA Typing"""

    #: Step name
    name = "hla_typing"
    produces = [DataSignature(DataType.TABULAR, frozenset({"hla"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    #: config_model_class
    config_model_class = HlaTypingConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA typing output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"txt": f"output/{lib}/out/{lib}.txt", "calls_json": f"output/{lib}/out/{lib}.json"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        match self.config.tool:
            case "optitype":
                selected = OptiTypeStepPart
            case "arcashla":
                selected = ArcasHlaStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes((LinkOutStepPart, selected))
        #: Mapping from library name to library object
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        from os.path import join

        name_pattern = "{library_name}"
        yield from self._yield_result_files(
            join("output", name_pattern, "out", name_pattern + "{ext}"), ext=EXT_VALUES
        )
        log_ext = [e + m for e in ("log", "conda_list.txt", "conda_info.txt") for m in ("", ".md5")]
        yield from self._yield_result_files(
            join("output", name_pattern, "log", name_pattern + ".{ext}"), ext=log_ext
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        tool = str(self.config.tool)
        supported = self.sub_steps[tool].supported_extraction_types
        df = self.build_library_dataframe()
        if df.empty:
            return
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "").lower()
            if extraction_type in supported:
                yield from expand(
                    tpl,
                    library_name=[library_name],
                    **kwargs,
                )
