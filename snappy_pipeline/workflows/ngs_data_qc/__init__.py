# -*- coding: utf-8 -*-
"""Implementation of the ``ngs_data_qc`` step

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_ngs_data_qc.rst

"""

import os
from itertools import chain

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import Namedlist, expand, touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
    get_ngs_library_folder_name,
)

#: Default configuration for the ngs_mapping schema
DEFAULT_CONFIG = r"""
# Default configuration ngs_mapping
step_config:
  ngs_data_qc:
    tools:
    - fastqc
"""


class FastQcReportStepPart(BaseStepPart):
    """(Raw) data QC using FastQC"""

    name = "fastqc"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_set_infos, self.parent.config_lookup_paths
        )

    def get_args(self, action):
        def args_function(wildcards):
            return {
                "num_threads": 1,
                "more_reads": Namedlist(
                    chain(
                        sorted(self._collect_reads(wildcards, wildcards.library_name, "")),
                        sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")),
                    )
                ),
            }

        assert action == "run", "Unsupported actions"
        return args_function

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            return "work/input_links/{library_name}/.done".format(**wildcards)

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files for the (raw) data QC steps"""
        assert action == "run"
        yield "fastqc_done", touch("work/{library_name}/report/fastqc/.done")

    def get_log_file(self, action):
        return "work/{library_name}/log/snakemake.fastqc.log"

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)

    def update_cluster_config(self, cluster_config):
        cluster_config["data_qc_fastqc_run"] = {
            "mem": int(3.75 * 1024 * 2),
            "time": "12:00",
            "ntasks": 2,
        }


class NgsDataQcWorkflow(BaseStep):
    """Perform NGS raw data QC"""

    name = "ngs_data_qc"
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
        )
        self.register_sub_step_classes((LinkInStep, LinkOutStepPart, FastQcReportStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS raw data QC workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        from os.path import join

        token = "{ngs_library.name}"
        # TODO: actually link out report files
        yield from self._yield_result_files(join("output", token, "report", "fastqc", ".done"))

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                # extraction_type = ngs_library.test_sample.extra_infos['extractionType']
                yield from expand(tpl, ngs_library=[ngs_library], **kwargs)
