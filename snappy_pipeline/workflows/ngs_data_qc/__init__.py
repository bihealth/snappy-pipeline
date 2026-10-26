# -*- coding: utf-8 -*-
"""Quality control of raw sequencing data

The step runs a QC tool on the reads or the alignments of each library.
The ``tool`` field selects the tool:

- ``fastqc``: runs FastQC on the FASTQ files of DNA and RNA libraries. It reads ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
- ``picard``: runs the Picard programs listed in ``picard.programs`` on the BAM files of DNA libraries. It reads ``depends_on.alignments``, for example an :ref:`step_ngs_mapping` task.
- ``ngs_chew``: writes a fingerprint of each DNA library with `ngs-chew <https://github.com/bihealth/ngs-chew>`_, its genotypes at selected common germline variants, from the BAM files of ``depends_on.alignments``.
  One job compares the fingerprints of all libraries of the task.

The result is one set of files per library.
The task processes all libraries, normals and tumors alike, unless ``library_selection`` narrows them.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.

With ``fastqc``:

- ``{library_name}/report/.done``, which marks the finished run. The FastQC HTML reports are in ``work/{library_name}/report/html/`` and are no result files.

With ``picard``, for each program in ``picard.programs``, with an ``.md5`` file:

- ``{library_name}/report/{library_name}.<program>.txt``
- For the programs that ``CollectMultipleMetrics`` runs, ``{library_name}/report/{library_name}.CollectMultipleMetrics.<first metrics file>.txt``, for example ``alignment_summary_metrics`` for ``CollectAlignmentSummaryMetrics``.

With ``ngs_chew``, with an ``.md5`` file:

- ``{library_name}/out/{library_name}.npz``, the fingerprint
- ``cohort/out/cohort.relatedness.txt``, ``cohort.ibs0.txt`` and ``cohort.mask.txt``, for each pair of libraries the relatedness, the number of sites where they share no allele (IBS0) and the number of sites both cover
- ``cohort/out/cohort.stats.txt``, for each library its heterozygous and homozygous sites, the variance of the allele fractions at its heterozygous sites (``var_het``) and its chrX and chrY fractions

Downstream tasks read the key ``done``, or ``npz`` with ``ngs_chew``.

Notes
-----

The Picard metrics are registered as Snakemake ``report()`` outputs.
With ``ngs_chew``, the report shows the relatedness of all pairs of libraries and the var(het) of each library.
Libraries of one donor have a relatedness close to 1 and those of unrelated donors close to 0, so a sample swap shows as a pair that breaks this pattern.
A library contaminated with DNA of another person has a higher var(het) than the others.
Fingerprints assign a sample to a patient more sensitively than HLA typing, which is matched between donor and recipient in stem-cell treatments.
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

#: The tables of ``ngs-chew compare`` and ``ngs-chew stats``
COMPARE_TABLES = ("relatedness", "ibs0", "mask", "stats")


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
        yield "baits", "work/{library_kit}/out/{library_kit}.baits.interval_list"
        yield "targets", "work/{library_kit}/out/{library_kit}.targets.interval_list"

    @dictify
    def _get_output_files_metrics(self):
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


class NgsChewStepPart(BaseStepPart):
    """Fingerprints with ngs-chew, compared across the libraries of the task"""

    name = "ngs_chew"
    actions = ("fingerprint", "compare")

    #: ``compare`` runs ngs-chew's comparison, statistics and plots, which take no
    #: ``extra_arguments``
    extra_arguments_action = "fingerprint"

    report_outputs = {
        "compare": {
            "relatedness_html": ReportOutput("ngs-chew", caption="ngs_chew_relatedness.rst"),
            "var_het_html": ReportOutput("ngs-chew", caption="ngs_chew_var_het.rst"),
        }
    }

    resource_usage = {"fingerprint": ResourceUsage(threads=1, runtime="4h", mem="2GB")}

    @dictify
    def _get_input_files_fingerprint(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield from self.parent.reference_inputs().items()

    def _get_output_files_fingerprint(self):
        npz = "work/{library_name}/out/{library_name}.npz"
        return {"npz": npz, "npz_md5": npz + ".md5"}

    def _get_params_fingerprint(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}  # get_params adds the extra_arguments

    def _get_log_file_fingerprint(self):
        return {"log": "work/{library_name}/log/{library_name}.fingerprint.log"}

    def _get_input_files_compare(self, wildcards):
        return {"npz": [f"work/{lib}/out/{lib}.npz" for lib in self.parent.dna_libraries()]}

    @dictify
    def _get_output_files_compare(self):
        for table in COMPARE_TABLES:
            yield table, f"work/cohort/out/cohort.{table}.txt"
            yield f"{table}_md5", f"work/cohort/out/cohort.{table}.txt.md5"
        yield "relatedness_html", "work/cohort/report/cohort.relatedness.html"
        yield "var_het_html", "work/cohort/report/cohort.var_het.html"

    def _get_log_file_compare(self):
        return {"log": "work/cohort/log/cohort.compare.log"}


class NgsDataQcWorkflow(BaseStep):
    """Perform NGS raw data QC"""

    name = "ngs_data_qc"
    config_model_class = NgsDataQcConfigModel
    produces = [DataSignature(DataType.QC)]

    #: Every library, not only the tumors of cancer sheets
    library_default = ""

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local NGS QC output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        if config.tool == "ngs_chew":
            return {"npz": f"work/{lib}/out/{lib}.npz"}
        return {"done": f"work/{lib}/report/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        tool_to_class = {
            "fastqc": FastQcReportStepPart,
            "picard": PicardStepPart,
            "ngs_chew": NgsChewStepPart,
        }
        self.register_sub_step_classes((tool_to_class[self.config.tool],))
        if self.config.tool == "picard" and "CollectHsMetrics" in self.config.picard.programs:
            self.check_library_kits("kits", self.dna_libraries())

    def dna_libraries(self) -> list[str]:
        """Return the DNA libraries that ``picard`` and ``ngs_chew`` run on."""
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
        if self.config.tool == "ngs_chew":
            yield from self._yield_result_files(
                tpl="output/{library_name}/out/{library_name}.{ext}",
                allowed_extraction_types=("dna",),
                ext=["npz", "npz.md5"],
            )
            for table in COMPARE_TABLES:
                yield f"output/cohort/out/cohort.{table}.txt"
                yield f"output/cohort/out/cohort.{table}.txt.md5"

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
