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
from collections import OrderedDict

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

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

#: Extensions of files to create as main payload
EXT_VALUES = (".txt", ".txt.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("txt", "txt_md5")

#: HLA typing tools
HLA_TYPERS = ("optitype", "arcashla")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = r"""
# Default configuration ngs_mapping
step_config:
  hla_typing:
    path_ngs_mapping: ../ngs_mapping
    path_link_in: ""
    tools:
      - optitype
    optitype:
      max_reads: 5000  # suggestion by OptiType author
      num_mapping_threads: 4
    arcashla:
      mapper: star
""".lstrip()


class OptiTypeStepPart(BaseStepPart):
    """HLA Typing using OptiType"""

    #: Step name
    name = "optitype"

    #: Class available actions
    actions = ("run",)

    #: Supported extraction types: DNA, RNA
    supported_extraction_types = ["dna", "rna"]

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
            preprocessed_path=self.config["path_link_in"],
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

    @staticmethod
    def get_log_file(action):
        """Return path to log file"""
        _ = action
        return "work/optitype.{library_name}/log/snakemake.hla_typing.log"

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
            return result

        assert action == "run", "Unsupported actions"
        return args_function

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        _ = library_name
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)

    def _get_seq_type(self, wildcards):
        """Return sequence type for the library name in wildcards"""
        library = self.parent.ngs_library_name_to_ngs_library[wildcards.library_name]
        return library.test_sample.extra_infos["extractionType"].lower()

    def get_resource_usage(self, action):
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

    #: Supported extraction type: RNA
    supported_extraction_types = ["rna"]

    def __init__(self, parent):
        super().__init__(parent)
        self.mapper = self.config["arcashla"]["mapper"]
        self.base_path_out = (
            "work/{mapper}.arcashla.{{library_name}}/out/{mapper}.arcashla.{{library_name}}{ext}"
        )
        self.extensions = EXT_VALUES

    def get_input_files(self, action):
        """Return input files"""

        @dictify
        def input_function(wildcards):
            yield "ref_done", "work/arcashla.prepare_reference/out/.done"
            tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
            yield "bam", self.parent.sub_workflows["ngs_mapping"](
                tpl.format(mapper=self.mapper, **wildcards)
            )

        assert action == "run"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext, mapper=self.config["arcashla"]["mapper"])

    def get_output_prefix(self):
        return "%s." % self.config["arcashla"]["mapper"]

    @staticmethod
    def get_log_file(action):
        """Return path to log file"""
        _ = action
        return "work/arcashla.{library_name}/log/snakemake.hla_typing.log"

    def get_resource_usage(self, action):
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
        super().__init__(*args, **kwargs)
        self.register_sub_step_classes(
            (OptiTypeStepPart, ArcasHlaStepPart, LinkInStep, LinkOutStepPart)
        )
        #: Mapping from library name to library object
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

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

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                for tool in self.config["tools"]:
                    supported = self.sub_steps[tool].supported_extraction_types
                    extraction_type = ngs_library.test_sample.extra_infos["extractionType"].lower()
                    if extraction_type in supported:
                        yield from expand(
                            tpl,
                            prefix=self.sub_steps[tool].get_output_prefix(),
                            hla_typer=[tool],
                            ngs_library=[ngs_library],
                            **kwargs,
                        )
