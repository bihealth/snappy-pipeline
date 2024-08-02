# -*- coding: utf-8 -*-
"""Implementation of the ``adapter_trimming`` step

The adapter_trimming step performs adapter & quality trimming of reads (DNA or RNA).
The tools are highly configurable, and provide feedback of the success of the operation.

==========
Step Input
==========

For each library defined in all sample sheets, the instances of this step will search for the input
files according to the configuration.  The found read files will be linked into
``work/input_links/{library_name}`` (status quo, not a output path, thus path not guaranteed
to be stable between minor versions).

The search paths can be overridden using the step configuration option ``path_link_in``.
``path_link_in`` is a general features that enables pre-processing steps, typically before mapping.

----------------------
Data Set Configuration
----------------------

Consider the following data set definition from the main configuration file.

.. code-block:: yaml

    data_sets:
      first_batch:
        file: 01_first_batch.tsv
        search_patterns:
          # Note that currently only "left" and "right" key known
          - {'left': '*/L???/*_R1.fastq.gz', 'right': '*/L???/*_R2.fastq.gz'}
        search_paths: ['../input/01_first_batch']

Here, the data set ``first_batch`` is defined.  The sample sheet file is named
``01_first_batch.tsv`` and looked for in the relative path to the configuration file.  The input
search will be start in the (one, but could be more than one) path ``../input/01_first_batch``
(relative to the directory containing the configuration file).  The sample sheet provides a
``folderName`` ``extraInfo`` entry for each NGS library.  This folder name is searched for (e.g.,
``P001-N1-DNA1-WES``).  Once such a folder is found, the patterns in the values of the dict
``search_patterns`` are used for locating the paths of the actual files.

Currently, the only supported keys in the ``search_patterns`` dict are ``"left"`` and ``"right""``
(the latter can be omitted when only searching for single-end reads).

Consider the following example:

::

  ../input/
    `-- 01_first_batch
        |-- P001-N1-DNA1-WES1
        |   `-- 42KF5AAXX
        |       `-- L001
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz.md5
        |           |-- P001-N1-DNA1-WES1_R2.fastq.gz
        |           `-- P001-N1-DNA1-WES1_R2.fastq.gz.md5
        [...]

Here, the folder ``01_first_batch`` will be searched for a directory named ``P001-N1-DNA1-WES``.
After finding, the relative paths ``42KF5AAXX/L001/P001-N1-DNA1-WES1_R1.fastq.gz`` and
``42KF5AAXX/L001/P001-N1-DNA1-WES1_R2.fastq.gz`` will be found and used for the left/right parts of
a paired read set.

------------------------------------------------------
Overriding data set confguration with ``path_link_in``
------------------------------------------------------

When the config option ``path_link_in`` is set, it takes precedence on the search paths defined in the
data set configuration.

The searching for input files will follow the same rules as defined in the data set configuration,
except that the base path for the search provided by one single path defined in the configuration of the
step.


Mixing Single-End and Paired-End Reads
======================================

By default, it is checked that for each ``search_pattern``, the same number of matching files
has to be found, otherwise directories are ignored.  The reason is to reduce the number of
possible errors when linking in files.  You can change this behaviour by specifying
``mixed_se_pe: True`` in the data set information.  Then, it will be allowed to have the matches
for the ``right`` entry to be empty.  However, you will need to consistently have either SE or
PE data for each library; it is allowed to mix SE and PE libraries within one project but not
to have PE and SE data for one library.

Note that mixing single-end and paired-end reads is not (yet) supported when overriding the data set
configuration by setting a value to the configuration option ``path_link_in``.


===========
Step Output
===========

Adapter trimming will be performed for all NGS libraries in all sample sheets.  For each combination
of tool library, a directory ``{tool}/{lib_name}-{lib_pk}/out`` will be created.
Therein, trimmed fastq files will be created.

The input structure and file names will be maintained on output.  For example, it might look as
follows for the example from above:

::

    output/
    +-- bbduk
    |   `-- out
    |       `-- P001-N1-DNA1-WES1
    |           |-- 42KF5AAXX
    |           |   `-- L001
    |           |       |-- P001-N1-DNA1-WES1_R1.fastq.gz
    |           |       |-- P001-N1-DNA1-WES1_R1.fastq.gz.md5
    |           |       |-- P001-N1-DNA1-WES1_R2.fastq.gz
    |           |       `-- P001-N1-DNA1-WES1_R2.fastq.gz.md5
    |           `-- .done
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

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    ResourceUsage,
    get_ngs_library_folder_name,
)

from .model import AdapterTrimming as AdapterTrimmingConfigModel

#: Adatper trimming tools
TRIMMERS = ("bbduk", "fastp")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = AdapterTrimmingConfigModel.default_config_yaml_string()


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
            preprocessed_path=self.config.path_link_in,
        )

    @dictify
    def get_input_files(self, action):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        yield "done", "work/input_links/{library_name}/.done"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        return (
            ("out_done", self.base_path_out.format(trimmer=self.name) + "/out/.done"),
            ("report_done", self.base_path_out.format(trimmer=self.name) + "/report/.done"),
            ("rejected_done", self.base_path_out.format(trimmer=self.name) + "/rejected/.done"),
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        _ = action
        prefix = "work/{trimmer}.{{library_name}}/log/{trimmer}.{{library_name}}".format(
            trimmer=self.__class__.name
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        yield (
            "done",
            "work/{trimmer}.{{library_name}}/log/.done".format(trimmer=self.__class__.name),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
            if self.config.path_link_in:
                folder_name = wildcards.library_name
            reads_left = self._collect_reads(wildcards, folder_name, "")
            reads_right = self._collect_reads(wildcards, folder_name, "right-")
            return {
                "library_name": wildcards.library_name,
                "input": {
                    "reads_left": {key: reads_left[key] for key in sorted(reads_left.keys())},
                    "reads_right": {key: reads_right[key] for key in sorted(reads_right.keys())},
                },
            }

        # Validate action
        self._validate_action(action)
        return args_function

    def _collect_reads(self, wildcards, folder_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        path_info = {}
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            input_path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            assert input_path not in path_info.keys()
            paths = {
                "relative_path": path_infix,
                "filename": filename,
            }
            path_info[input_path] = paths
        return path_info


class BbdukStepPart(AdapterTrimmingStepPart):
    """bbduk adapter trimming"""

    name = "bbduk"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.bbduk.num_threads,
            time="12:00:00",  # 40 hours
            memory="24000M",
        )


class FastpStepPart(AdapterTrimmingStepPart):
    """fastp adapter trimming"""

    #: Step name
    name = "fastp"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.fastp.num_threads,
            time="12:00:00",  # 60 hours
            memory="24000M",
        )


class LinkOutFastqStepPart(BaseStepPart):
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
            """Helper wrapper function"""
            return expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)

        self._validate_action(action)
        return input_function

    def get_output_files(self, action):
        self._validate_action(action)
        return expand(self.base_path_out, sub_dir=self.sub_dirs)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out postprocessed (or not) files"""
        # Validate action
        self._validate_action(action)
        ins = expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)
        outs = [s.format(**wildcards) for s in expand(self.base_path_out, sub_dir=self.sub_dirs)]
        assert len(ins) == len(outs)

        cmd = "din_=$(dirname {in_}) ; dout=$(dirname {out})"
        cmd = cmd + " ; fns=$(find $din_ -type f -printf '%P\\n')"
        cmd = (
            cmd
            + " ; for fn in $fns ; do"
            + "     if [[ ! -L $din_/$fn ]] ; then"
            + "       mkdir -p $(dirname $dout/$fn) ; ln -sr $din_/$fn $dout/$fn"
            + "   ; fi"
            + " ; done"
        )
        return "\n".join((cmd.format(in_=in_, out=out) for in_, out in zip(ins, outs)))

    def _validate_action(self, action):
        assert action == "run"


class AdapterTrimmingWorkflow(BaseStep):
    """Perform adapter & quality-based trimming"""

    #: Step name
    name = "adapter_trimming"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, config_model_class=AdapterTrimmingConfigModel)
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
                for tool in self.config.tools:
                    for tpl in tpls:
                        yield tpl.format(trimmer=tool, ngs_library_name=ngs_library.name)
