# -*- coding: utf-8 -*-
"""Implementation of the ``adapter_trimming`` step

The adapter_trimming step performs adapter & quality trimming of reads (DNA or RNA).
The tools are highly configurable, and provide feedback of the success of the operation.

==========
Step Input
==========

Adapter trimming starts at the raw RNA-seq reads.  Thus, the input is very similar to one of
:ref:`ngs_mapping step <step_ngs_mapping>`.

See :ref:`ngs_mapping_step_input` for more information.

===========
Step Output
===========

Adapter trimming will be performed for all NGS libraries in all sample sheets.  For each combination
of tool library, a directory ``{tool}.{lib_name}-{lib_pk}/out`` will be created.
Therein, trimmed fastq files will be created, together with fastq files for reads that might be
filtered out by the process.

For example, it might look as follows for the example from above:

::

    output/
    +-- bbduk.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- Sample001_Lane_001_R1.fastq.gz
    |       |-- Sample001_Lane_001_R2.fastq.gz
    |       |-- Sample001_Lane_002_R1.fastq.gz
    |       |-- Sample001_Lane_002_R2.fastq.gz
    |       |-- Sample001_Lane_001_R1.fastq.gz.md5
    |       |-- Sample001_Lane_001_R2.fastq.gz.md5
    |       |-- Sample001_Lane_002_R1.fastq.gz.md5
    |       |-- Sample001_Lane_002_R2.fastq.gz.md5
    |       `-- .done
    [...]

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_adapter_trimming.rst

================================
Available Adapter Trimming Tools
================================

The following adpter trimming tools are currently available

- ``"bbduk"``
- ``"fastp"``

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

#: Adatper trimming tools
TRIMMERS = ("bbduk", "fastp")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = r"""
# Default configuration ngs_mapping
step_config:
  adapter_trimming:
    path_link_in: ""
    tools:
      - bbduk
    bbduk:
      adapter_sequences: REQUIRED # REQUIRED
      num_threads: 8
      parameters:  # Taken from https://github.com/ewels/MultiQC/issues/1146#issuecomment-607980076
      - "minlength=35"
      - "trimq=25"
      - "ktrim=r"
      - "qtrim=rl"
      - "k=18"
      - "mink=11"
      - "hdist=1"
      - "trimpolyg=8" # Remove long polyG
    fastp:
      num_threads: 4
      parameters:
      - "--trim_poly_g"
      - "--poly_g_min_len 8"
""".lstrip()


class AdapterTrimmingStepPart(BaseStepPart):
    """Adapter trimming common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{trimmer}.{{library_name}}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

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
        return (
            ("out_done", self.base_path_out.format(trimmer=self.name) + "/out/.done"),
            ("report_done", self.base_path_out.format(trimmer=self.name) + "/report/.done"),
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = "work/{trimmer}.{{library_name}}/log/{trimmer}.{{library_name}}".format(
            trimmer=self.__class__.name
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        yield "done", "work/{trimmer}.{{library_name}}/log/.done".format(
            trimmer=self.__class__.name
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


class BbdukStepPart(AdapterTrimmingStepPart):
    """bbduk adapter trimming"""

    name = "bbduk"

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
            threads=self.config["bbduk"]["num_threads"],
            time="12:00:00",  # 40 hours
            memory="24000M",
        )


class FastpStepPart(AdapterTrimmingStepPart):
    """fastp adapter trimming"""

    #: Step name
    name = "fastp"

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
            threads=self.config["fastp"]["num_threads"],
            time="12:00:00",  # 60 hours
            memory="24000M",
        )


class LinkOutFastqStepPart(LinkOutStepPart):
    """Link out the trimming results (all fastqs)"""

    name = "link_out_fastq"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/{wildcards.trimmer}.{wildcards.library_name}/{{sub_dir}}/.done"
        self.base_path_out = "output/{{trimmer}}/{{library_name}}/{sub_dir}/.done"
        self.sub_dirs = ["log", "report", "out"]

    def get_input_files(self, action):
        """Return required input files"""

        def input_function(wildcards):
            """Helper rapper function"""
            return expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)

        assert action == "run", "Unsupported action"
        return input_function

    def get_output_files(self, action):
        """Return output files that are generated by snappy-gatk_post_bam"""
        assert action == "run", "Unsupported action"
        return expand(self.base_path_out, sub_dir=self.sub_dirs)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out postprocessed (or not) files"""
        assert action == "run", "Unsupported action"
        ins = expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)
        outs = [s.format(**wildcards) for s in expand(self.base_path_out, sub_dir=self.sub_dirs)]
        assert len(ins) == len(outs)

        cmd = "din_=$(dirname {in_}) ; dout=$(dirname {out})"
        cmd = cmd + " ; fns=$(ls $din_ | grep -Ev '^(rejected|unpaired|failed)_')"
        cmd = (
            cmd
            + " ; for fn in $fns ; do if [[ ! -L $din_/$fn ]] ; then ln -sr $din_/$fn $dout/$fn ; fi ; done"
        )
        cmd = cmd + " ; test -L {in_} || ln -sr {in_} {out}"
        return "\n".join((cmd.format(in_=in_, out=out) for in_, out in zip(ins, outs)))


class AdapterTrimmingWorkflow(BaseStep):
    """Perform adapter & quality-based trimming"""

    #: Step name
    name = "adapter_trimming"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.register_sub_step_classes(
            (BbdukStepPart, FastpStepPart, LinkInStep, LinkOutFastqStepPart)
        )
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    @listify
    def get_result_files(self):
        """Return list of fixed name result files for the adapter trimming workflow"""
        tpls = (
            "output/{trimmer}/{ngs_library_name}/out/.done",
            "output/{trimmer}/{ngs_library_name}/report/.done",
            "output/{trimmer}/{ngs_library_name}/log/.done",
        )
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                for tool in self.config["tools"]:
                    for tpl in tpls:
                        yield tpl.format(trimmer=tool, ngs_library_name=ngs_library.name)
