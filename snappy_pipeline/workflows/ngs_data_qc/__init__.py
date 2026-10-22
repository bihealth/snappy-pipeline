# -*- coding: utf-8 -*-
"""Quality control of raw sequencing data

The step runs a QC tool on the reads or the alignments of each library.
The ``tool`` field selects the tool:

- ``fastqc``: runs FastQC on the FASTQ files of DNA and RNA libraries. It reads ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
- ``picard``: runs the Picard programs listed in ``picard.programs`` on the BAM files of DNA libraries. It reads ``depends_on.alignments``, for example an :ref:`step_ngs_mapping` task.

The result is one set of files per library.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.

With ``fastqc``:

- ``{library_name}/report/.done``, which marks the finished run. The FastQC HTML reports are in ``work/{library_name}/report/html/`` and are not published.

With ``picard``, for each program in ``picard.programs``, with an ``.md5`` file:

- ``{library_name}/report/{library_name}.<program>.txt``
- For the programs that ``CollectMultipleMetrics`` runs, ``{library_name}/report/{library_name}.CollectMultipleMetrics.<first metrics file>.txt``, for example ``alignment_summary_metrics`` for ``CollectAlignmentSummaryMetrics``.

Downstream tasks read the key ``done``.

Notes
-----

The Picard metrics are registered as Snakemake ``report()`` outputs.
``CollectHsMetrics`` needs ``depends_on.kits``, an ``external_data`` task with ``kits``, and a kit for each DNA library, so it runs on targeted libraries only.
The step converts the baits and targets BED files of each kit to interval lists in one preparation job per kit, ``work/{library_kit}/out/{library_kit}.{baits,targets}.interval_list``.
The kit name labels the ``BAIT_SET`` column of the metrics.
"""

from typing import Any

from snakemake.io import directory, expand, touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import NgsDataQc as NgsDataQcConfigModel


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
    def _get_output_files_run(self):
        yield "fastqc_done", touch("work/{library_name}/report/.done")
        # One report per FASTQ file, so a directory with an index of them
        yield "html", directory("work/{library_name}/report/html")

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}


class PicardStepPart(BaseStepPart):
    """Collect Picard metrics"""

    name = "picard"
    actions = ("prepare", "metrics")

    #: ``prepare`` converts the BED files of the kit to interval lists
    extra_arguments_action = "metrics"
    intermediate_outputs = {"prepare": ("baits", "targets")}

    @property
    def report_outputs(self):
        if self.config.tool != self.name:
            return {}
        programs = self.config.picard.programs
        metrics = {pgm: ReportOutput(f"Picard {pgm}", caption="picard.rst") for pgm in programs}
        return {"metrics": metrics}

    @dictify
    def _get_input_files_prepare(self, wildcards):
        yield from self.parent.reference_inputs().items()
        kit = self.parent.get_upstream_paths("kits", library_kit=wildcards.library_kit)
        yield "baits", kit.baits
        yield "targets", kit.targets

    @dictify
    def _get_input_files_metrics(self, wildcards):
        yield from self.parent.reference_inputs().items()
        if self.config.depends_on.dbsnp:
            yield from self.parent.vcf_inputs("dbsnp").items()
        if "CollectHsMetrics" in self.config.picard.programs:
            kit = self.parent.library_kits[wildcards.library_name]
            yield "baits", f"work/{kit}/out/{kit}.baits.interval_list"
            yield "targets", f"work/{kit}/out/{kit}.targets.interval_list"
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    @dictify
    def _get_output_files_prepare(self):
        if self.name != self.config.tool:
            return
        yield "baits", "work/{library_kit}/out/{library_kit}.baits.interval_list"
        yield "targets", "work/{library_kit}/out/{library_kit}.targets.interval_list"

    @dictify
    def _get_output_files_metrics(self):
        if self.name != self.config.tool:
            return
        base_out = "work/{library_name}/report/{library_name}."
        for pgm in self.config.picard.programs:
            if pgm in MULTIPLE_METRICS.keys():
                first = MULTIPLE_METRICS[pgm][0]
                yield pgm, base_out + f"CollectMultipleMetrics.{first}.txt"
                yield pgm + "_md5", base_out + f"CollectMultipleMetrics.{first}.txt.md5"
            else:
                yield pgm, base_out + pgm + ".txt"
                yield pgm + "_md5", base_out + pgm + ".txt.md5"

    def _get_log_file_prepare(self):
        return {"log": "work/{library_kit}/log/{library_kit}.prepare.log"}

    def _get_log_file_metrics(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_metrics(self, wildcards: Wildcards) -> dict[str, Any]:
        params = {
            "prefix": f"{wildcards.library_name}.",
            "programs": self.config.picard.programs,
        }
        if "CollectHsMetrics" in self.config.picard.programs:
            params["bait_name"] = self.parent.library_kits[wildcards.library_name]
        return params

    resource_usage = {"metrics": ResourceUsage(threads=1, runtime="24h", mem="64GB")}


class NgsDataQcWorkflow(BaseStep):
    """Perform NGS raw data QC"""

    name = "ngs_data_qc"
    config_model_class = NgsDataQcConfigModel
    produces = [DataSignature(DataType.QC)]

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local NGS QC output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"work/{lib}/report/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((FastQcReportStepPart, PicardStepPart))
        if self.config.tool == "picard" and "CollectHsMetrics" in self.config.picard.programs:
            self.check_library_kits("kits", self._picard_libraries())

    def _picard_libraries(self) -> list[str]:
        """Return the DNA libraries that ``picard`` runs on."""
        df = self.build_library_dataframe()
        dna = df[df["extraction_type"] == "dna"]["library_name"]
        return [name for name in self.output_entities if name in set(dna)]

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
