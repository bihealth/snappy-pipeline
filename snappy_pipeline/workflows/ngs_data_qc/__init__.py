# -*- coding: utf-8 -*-
"""Implementation of the ``ngs_data_qc`` step

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_ngs_data_qc.rst

"""

from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import directory, expand, touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import NgsDataQc as NgsDataQcConfigModel

#: Default configuration for the ngs_mapping schema

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

    default_resource_usage = ResourceUsage(threads=1, mem="4GB", runtime="4h")

    report_outputs = {
        "run": {"html": ReportOutput("FastQC", caption="fastqc.rst", htmlindex="index.html")}
    }

    def _get_params_run(self, wildcards):
        return {"num_threads": 1} | reads_params(self.parent, wildcards.library_name)

    def _get_input_files_run(self, wildcards):
        return reads_input_files(self.parent, wildcards.library_name)

    @dictify
    def get_output_files(self, action):
        """Return output files for the (raw) data QC steps"""
        # Validate action
        self._validate_action(action)
        yield "fastqc_done", touch("work/{library_name}/report/.done")
        # One report per FASTQ file, so a directory with an index of them
        yield "html", directory("work/{library_name}/report/html")

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        prefix = "work/{library_name}/log/{library_name}"
        yield "log", prefix + ".log"
        yield "conda_info", prefix + ".conda_info.txt"
        yield "conda_list", prefix + ".conda_list.txt"


class PicardStepPart(BaseStepPart):
    """Collect Picard metrics"""

    name = "picard"
    actions = ("prepare", "metrics")

    @property
    def report_outputs(self):
        if self.config.tool != self.name:
            return {}
        programs = self.config.picard.programs
        metrics = {pgm: ReportOutput(f"Picard {pgm}", caption="picard.rst") for pgm in programs}
        return {"metrics": metrics}

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_prepare(self, wildcards):
        yield "reference", self.parent.get_upstream_paths("reference").fasta
        if self.config.picard.path_to_baits:
            yield "baits", self.config.picard.path_to_baits
        if self.config.picard.path_to_targets:  # same as the baits when missing
            yield "targets", self.config.picard.path_to_targets

    @dictify
    def _get_input_files_metrics(self, wildcards):
        yield "reference", self.parent.get_upstream_paths("reference").fasta
        if self.config.depends_on.dbsnp:
            yield "dbsnp", self.parent.get_upstream_paths("dbsnp").vcf
        if "CollectHsMetrics" in self.config.picard.programs:
            yield "baits", "work/static_data/out/baits.interval_list"
            yield "targets", "work/static_data/out/targets.interval_list"
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    @dictify
    def get_output_files(self, action):
        if self.name != self.config.tool:
            return {}
        if action == "prepare":
            yield "baits", "work/static_data/out/baits.interval_list"
            yield "targets", "work/static_data/out/targets.interval_list"
        elif action == "metrics":
            base_out = "work/{library_name}/report/{library_name}."
            for pgm in self.config.picard.programs:
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
            prefix = "work/static_data/log/prepare"
        elif action == "metrics":
            prefix = "work/{library_name}/log/{library_name}"
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

    def _get_params_metrics(self, wildcards: Wildcards) -> dict[str, Any]:
        params = {
            "prefix": f"{wildcards.library_name}.",
            "programs": self.config.picard.programs,
        }
        if self.config.picard.bait_name:
            params["bait_name"] = self.config.picard.bait_name
        return params

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action == "prepare":
            return super().get_resource_usage(action, **kwargs)
        elif action == "metrics":
            return ResourceUsage(threads=1, runtime="24h", mem="64GB")
        else:
            actions_str = ", ".join(self.actions)
            raise UnsupportedActionException(
                f"Action '{action}' is not supported. Valid options: {actions_str}"
            )


class NgsDataQcWorkflow(BaseStep):
    """Perform NGS raw data QC"""

    name = "ngs_data_qc"
    config_model_class = NgsDataQcConfigModel
    produces = [DataSignature(DataType.QC)]
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local NGS QC output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"output/{lib}/report/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((LinkOutStepPart, FastQcReportStepPart, PicardStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS raw data QC workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if self.config.tool == "fastqc":
            yield from self._yield_result_files(
                tpl="output/{library_name}/report/.done",
                allowed_extraction_types=(
                    "dna",
                    "rna",
                ),
            )
        if self.config.tool == "picard":
            tpl = "output/{library_name}/report/{library_name}.{ext}"
            exts = []
            for pgm in self.config.picard.programs:
                if pgm in MULTIPLE_METRICS.keys():
                    first = MULTIPLE_METRICS[pgm][0]
                    exts.append(f"CollectMultipleMetrics.{first}.txt")
                    exts.append(f"CollectMultipleMetrics.{first}.txt.md5")
                else:
                    exts.append(pgm + ".txt")
                    exts.append(pgm + ".txt.md5")
            yield from self._yield_result_files(
                tpl=tpl,
                allowed_extraction_types=("dna",),
                ext=exts,
            )

    def _yield_result_files(self, tpl, allowed_extraction_types, **kwargs):
        """Build output paths from path template and extension list"""
        df = self.build_library_dataframe()
        if df.empty:
            return
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type in allowed_extraction_types:
                yield from expand(tpl, library_name=[library_name], **kwargs)
