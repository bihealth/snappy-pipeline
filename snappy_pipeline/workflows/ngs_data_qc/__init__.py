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

#: Default configuration for the ngs_mapping schema
DEFAULT_CONFIG = r"""
# Default configuration ngs_mapping
step_config:
  ngs_data_qc:
    path_link_in: ""  # OPTIONAL Override data set configuration search paths for FASTQ files
    tools: [fastqc, picard]  # REQUIRED - available: 'fastqc' & 'picard' (for QC on bam files)
    picard:
      path_ngs_mapping: ../ngs_mapping  # REQUIRED
      path_to_baits: ""                 # Required when CollectHsMetrics is among the programs
      path_to_targets: ""               # When missing, same as baits
      bait_name: ""                     # Exon enrichment kit name (optional)
      programs: []  # Available metrics:
                    # * Generic metrics [* grouped into CollectMultipleMetrics]
      #                 - CollectAlignmentSummaryMetrics      *
      #                 - CollectBaseDistributionByCycle      *
      #                 - CollectGcBiasMetrics                *
      #                 - CollectInsertSizeMetrics            *
      #                 - CollectJumpingLibraryMetrics
      #                 - CollectOxoGMetrics
      #                 - CollectQualityYieldMetrics          *
      #                 - CollectSequencingArtifactMetrics    *
      #                 - EstimateLibraryComplexity
      #                 - MeanQualityByCycle                  *
      #                 - QualityScoreDistribution            *
                    # * WGS-specific metrics
      #                 - CollectRawWgsMetrics
      #                 - CollectWgsMetrics
      #                 - CollectWgsMetricsWithNonZeroCoverage
                    # * Other assay-specific metrics
      #                 - CollectHsMetrics                    Whole Exome Sequencing
      #                 - CollectTargetedPcrMetrics           Panel sequencing
      #                 - CollectRnaSeqMetrics                mRNA sequencing, not implemented yet
      #                 - CollectRbsMetrics                   bi-sulfite sequencing, not implemented yet
"""

MULTIPLE_METRICS = {
    "CollectAlignmentSummaryMetrics": ["alignment_summary_metrics"],
    "CollectBaseDistributionByCycle": ["base_distribution_by_cycle_metrics"],
    "CollectGcBiasMetrics": ["gc_bias.summary_metrics", "gc_bias.detail_metrics"],
    "CollectInsertSizeMetrics": ["insert_size_metrics"],
    "CollectQualityYieldMetrics": ["quality_yield_metrics"],
    "CollectSequencingArtifactMetrics": [
        "pre_adapter_detail_metrics",
        "pre_adapter_summary_metrics",
        "bait_bias_summary_metrics",
        "bait_bias_detail_metrics",
    ],
    "MeanQualityByCycle": ["quality_by_cycle_metrics"],
    "QualityScoreDistribution": ["quality_distribution_metrics"],
}
ADDITIONAL_METRICS = (
    "CollectJumpingLibraryMetrics",
    "CollectOxoGMetrics",
    "EstimateLibraryComplexity",
)
WGS_METRICS = (
    "CollectRawWgsMetrics",
    "CollectWgsMetrics",
    "CollectWgsMetricsWithNonZeroCoverage",
)
WES_METRICS = ("CollectHsMetrics",)
PANEL_METRICS = ("CollectTargetedPcrMetrics",)
RNA_METRICS = ("CollectRnaSeqMetrics",)
BISULFITE_METRICS = ("CollectRbsMetrics",)


class FastQcReportStepPart(BaseStepPart):
    """(Raw) data QC using FastQC"""

    #: Step name
    name = "fastqc"

    #: Class available actions
    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, memory="4G", time="03:59:59")

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

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

        return args_function

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            """Helper wrapper function"""
            return "work/input_links/{library_name}/.done".format(**wildcards)

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files for the (raw) data QC steps"""
        # Validate action
        self._validate_action(action)
        yield "fastqc_done", touch("work/{library_name}/report/fastqc/.done")

    @staticmethod
    def get_log_file(action):
        _ = action
        return "work/{library_name}/log/snakemake.fastqc.log"

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config["path_link_in"]:
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)


class PicardStepPart(BaseStepPart):
    """Collect Picard metrics"""

    name = "picard"
    actions = ("prepare", "metrics")

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        self._validate_action(action)
        if action == "prepare":
            raise UnsupportedActionException(
                'Action "prepare" input files must be defined in config'
            )

        return self._get_input_files_metrics

    @dictify
    def _get_input_files_metrics(self, wildcards):
        if "CollectHsMetrics" in self.config["picard"]["programs"]:
            yield "baits", "work/static_data/picard/out/baits.interval_list"
            yield "targets", "work/static_data/picard/out/targets.interval_list"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def get_output_files(self, action):
        if action == "prepare":
            yield "baits", "work/static_data/picard/out/baits.interval_list"
            yield "targets", "work/static_data/picard/out/targets.interval_list"
        elif action == "metrics":
            base_out = "work/{mapper}.{library_name}/report/picard/{mapper}.{library_name}."
            for pgm in self.config["picard"]["programs"]:
                if pgm in MULTIPLE_METRICS.keys():
                    first = MULTIPLE_METRICS[pgm][0]
                    yield pgm, base_out + f"CollectMultipleMetrics.{first}.txt"
                    yield pgm + "_md5", base_out + f"CollectMultipleMetrics.{first}.txt.md5"
                else:
                    yield pgm, base_out + pgm + ".txt"
                    yield pgm + "_md5", base_out + pgm + ".txt.md5"
        else:
            actions_str = ", ".join(self.actions)
            raise UnsupportedActionException(
                f"Action '{action}' is not supported. Valid options: {actions_str}"
            )

    @dictify
    def get_log_file(self, action):
        if action == "prepare":
            prefix = "work/static_data/picard/log/prepare"
        elif action == "metrics":
            prefix = "work/{mapper}.{library_name}/log/picard/{mapper}.{library_name}"
        else:
            actions_str = ", ".join(self.actions)
            raise UnsupportedActionException(
                f"Action '{action}' is not supported. Valid options: {actions_str}"
            )

        key_ext = (
            ("wrapper", ".wrapper.py"),
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_params(self, action):
        self._validate_action(action)

        return self._get_params

    @dictify
    def _get_params(self, wildcards):
        return {"prefix": f"{wildcards.mapper}.{wildcards.library_name}"}

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action == "prepare":
            return super().get_resource_usage(action)
        elif action == "metrics":
            return ResourceUsage(threads=1, time="24:00:00", memory="24G")
        else:
            actions_str = ", ".join(self.actions)
            raise UnsupportedActionException(
                f"Action '{action}' is not supported. Valid options: {actions_str}"
            )


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

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(workflow, config, config_lookup_paths, config_paths, workdir)
        self.register_sub_step_classes(
            (LinkInStep, LinkOutStepPart, FastQcReportStepPart, PicardStepPart)
        )
        if "picard" in self.config["tools"]:
            self.register_sub_workflow("ngs_mapping", self.config["picard"]["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS raw data QC workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if "fastqc" in self.config["tools"]:
            yield from self._yield_result_files(
                tpl="output/{ngs_library.name}/report/fastqc/.done",
                allowed_extraction_types=(
                    "DNA",
                    "RNA",
                ),
            )
        if "picard" in self.config["tools"]:
            tpl = (
                "output/{mapper}.{ngs_library.name}/report/picard/{mapper}.{ngs_library.name}.{ext}"
            )
            exts = []
            for pgm in self.config["picard"]["programs"]:
                if pgm in MULTIPLE_METRICS.keys():
                    first = MULTIPLE_METRICS[pgm][0]
                    exts.append(f"CollectMultipleMetrics.{first}.txt")
                    exts.append(f"CollectMultipleMetrics.{first}.txt.md5")
                else:
                    exts.append(pgm + ".txt")
                    exts.append(pgm + ".txt.md5")
            yield from self._yield_result_files(
                tpl=tpl,
                allowed_extraction_types=("DNA",),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                ext=exts,
            )

    def _yield_result_files(self, tpl, allowed_extraction_types, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                extraction_type = ngs_library.test_sample.extra_infos["extractionType"]
                if extraction_type in allowed_extraction_types:
                    yield from expand(tpl, ngs_library=[ngs_library], **kwargs)

    def check_config(self):
        if "picard" in self.config["tools"]:
            self.ensure_w_config(
                ("step_config", "ngs_data_qc", "picard", "path_ngs_mapping"),
                "Path to ngs_mapping not configured but required for picard",
            )
            programs = self.config["picard"]["programs"]
            assert len(programs) > 0, "No selected programs for collecting metrics"
            assert all(
                pgm in MULTIPLE_METRICS.keys()
                or pgm in ADDITIONAL_METRICS
                or pgm in WES_METRICS
                or pgm in WGS_METRICS
                for pgm in programs
            ), "Some requested metrics programs are not implemented"
            if "CollectHsMetrics" in programs:
                assert self.config["picard"][
                    "path_to_baits"
                ], "Path to baits must be specified when using CollectHsMetrics"
