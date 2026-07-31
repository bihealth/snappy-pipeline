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
from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand, touch
from snakemake.iocontainers import Namedlist, Wildcards

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStepPart,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

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

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.parent.get_preprocessed_path(),
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
        task_prefix = self.parent.task_path_prefix()
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.parent.get_preprocessed_path():
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            path = task_prefix + path
            yield path


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
        if "CollectHsMetrics" in self.config.picard.programs:
            yield "baits", "work/static_data/picard/out/baits.interval_list"
            yield "targets", "work/static_data/picard/out/targets.interval_list"
        infix = f"{wildcards.library_name}"
        yield (
            "bam",
            self.parent.upstream("ngs_mapping")(f"output/{infix}/out/{infix}.bam"),
        )

    @dictify
    def get_output_files(self, action):
        if self.name != self.config.tool:
            return {}
        if action == "prepare":
            yield "baits", "work/static_data/picard/out/baits.interval_list"
            yield "targets", "work/static_data/picard/out/targets.interval_list"
        elif action == "metrics":
            base_out = "work/{library_name}/report/picard/{library_name}."
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
            prefix = "work/static_data/picard/log/prepare"
        elif action == "metrics":
            prefix = "work/{library_name}/log/picard/{library_name}"
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

    def get_args(self, action):
        self._validate_action(action)

        return getattr(self, f"_get_args_{action}")

    def _get_args_prepare(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "path_to_baits": self.config.picard.path_to_baits,
            "path_to_targets": self.config.picard.path_to_targets,
        }

    def _get_args_metrics(self, wildcards: Wildcards) -> dict[str, Any]:
        params = {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "prefix": f"{wildcards.library_name}.",
            "programs": self.config.picard.programs,
        }
        if self.config.picard.bait_name:
            params["bait_name"] = self.config.picard.bait_name
        if (
            getattr(self.parent.w_config.static_data_config, "dbsnp", {"path": ""})
            and getattr(self.parent.w_config.static_data_config.dbsnp, "path", "")
            and self.parent.w_config.static_data_config.dbsnp.path
        ):
            params["dbsnp"] = self.parent.w_config.static_data_config.dbsnp.path
        else:
            params["dbsnp"] = ""
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
    consumes = {DataSignature(DataType.RAW): True, DataSignature(DataType.ALIGNMENTS): False}
    produces = [DataSignature(DataType.QC)]
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local NGS QC output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"output/{lib}/report/fastqc/.done"}

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            task_name=task_name,
            **kwargs,
        )
        self.register_sub_step_classes(
            (LinkInStepPart, LinkOutStepPart, FastQcReportStepPart, PicardStepPart)
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS raw data QC workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if self.config.tool == "fastqc":
            yield from self._yield_result_files(
                tpl="output/{library_name}/report/fastqc/.done",
                allowed_extraction_types=(
                    "DNA",
                    "RNA",
                ),
            )
        if self.config.tool == "picard":
            tpl = "output/{library_name}/report/picard/{library_name}.{ext}"
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
                allowed_extraction_types=("DNA",),
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
