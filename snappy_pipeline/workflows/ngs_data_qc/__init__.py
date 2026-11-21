# -*- coding: utf-8 -*-
"""Quality control of raw sequencing data

The step runs a QC tool on the reads or the alignments of each library.
The ``tool`` field selects the tool:

- ``fastqc``: runs FastQC on the FASTQ files of DNA and RNA libraries. It reads ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
- ``picard``: runs the Picard programs listed in ``picard.programs`` on the BAM files of DNA libraries. It reads ``depends_on.alignments``, for example an :ref:`step_ngs_mapping` task.
- ``ngs_chew``: writes a fingerprint of each DNA library with `ngs-chew <https://github.com/bihealth/ngs-chew>`_, its genotypes at selected common germline variants, from the BAM files of ``depends_on.alignments``.
  One job compares the fingerprints of all libraries of the task.
- ``alfred``: reports the coverage of the targets of each DNA library's kit with ``alfred qc``, for the libraries that have a kit (``WES``, ``Panel_seq``); it reads the ``targets`` of the kit from ``depends_on.kits``.
  ``varfish_export`` reads the report through its ``depends_on.qc``.

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

With ``alfred``, with an ``.md5`` file:

- ``{library_name}/out/{library_name}.target_coverage.json.gz``

With ``picard``, ``alfred`` and ``ngs_chew``, one job per task summarizes the files of all libraries in ``qc_summary/out/qc_summary.tsv``, one row per library (the QC summary of the notes).

Downstream tasks read the key ``done``, ``npz`` with ``ngs_chew`` or ``alfred_qc`` with ``alfred``.

Notes
-----

The Picard metrics are registered as Snakemake ``report()`` outputs.
With ``picard``, ``alfred`` and ``ngs_chew``, the QC summary job writes ``qc_summary/out/qc_summary.tsv`` and a site of tables and charts below ``work/qc_summary/report/html/``, which ``snappy report`` shows.
For ``picard``, the site has one table per metrics class and a chart of each histogram (target coverage, insert size, base quality by cycle and quality distribution); the summary table joins the metrics of each library into one row, the columns prefixed with the metrics class.
For ``alfred``, the summary table lists the mapping and coverage statistics of each library with the fraction of targets covered at least 10, 20, 30, 50 and 100 times, and the charts show the targets above a coverage, the on-target rate, the coverage and insert size distributions, the base quality by position and the GC content.
With ``ngs_chew``, the report shows the relatedness of all pairs of libraries and the var(het) of each library, and the fingerprint checks: one row per library with the sample sheet's sex against the sex predicted from the chrX genotypes (the sheet sex whose cohort median ratio is nearer, so the cohort needs both sexes), the sheet's father and mother against the relatedness and IBS0 rate of the pair, whether every pair of the pedigree has the relatedness the sheet implies (kinship over the pedigree, unsequenced donors included) and the libraries of other pedigrees that are related at all; a scatter chart of the libraries by chrX ratio and chrY reads colored by the sheet's sex, a scatter chart of the pairs (sheet against fingerprints) and the table of the pairs.
Libraries of one donor have a relatedness close to 1 and those of unrelated donors close to 0, so a sample swap shows as a pair that breaks this pattern.
A library contaminated with DNA of another person has a higher var(het) than the others.
Fingerprints assign a sample to a patient more sensitively than HLA typing, which is matched between donor and recipient in stem-cell treatments.
``CollectHsMetrics`` needs ``depends_on.kits``, an ``external_data`` task with ``kits``, and a kit for each DNA library, so it runs on targeted libraries only.
The step converts the baits and targets BED files of each kit to interval lists in one preparation job per kit, ``work/{library_kit}/out/{library_kit}.{baits,targets}.interval_list``.
The kit name labels the ``BAIT_SET`` column of the metrics.

Enrichment kits
~~~~~~~~~~~~~~~

``alfred`` and Picard's ``CollectHsMetrics`` read the region files of each library's enrichment kit.
The files come from an ``external_data`` task in ``kits`` mode, named in ``depends_on.kits``.
A library uses the kit named exactly like its ``libraryKit`` in the sample sheet:

::

    [Custom Fields]
    key     annotatedEntity docs    type    minimum maximum unit    choices pattern
    libraryKit      ngsLibrary      exome enrichment kit    string  0       0       0       0       0

    [Data]
    patientName     sampleName      isTumor extractionType  libraryType     folderName      libraryKit
    case001 N1      N       DNA     WES     case001-N1-DNA1-WES1    V5_noUTR
    case002 N1      N       DNA     WES     case002-N1-DNA1-WES1    V5_UTR
    [...]

Only targeted libraries (library type ``WES`` or ``Panel_seq``) have a kit; ``alfred`` skips all others.
A ``libraryKit`` that the kits task lacks is an error.

.. code-block:: yaml

    - step: external_data
      name: kits
      config:
        produces: {type: regions, tags: [kits]}
        kits:
          V5_noUTR:
            targets: resources/v5.bed
            baits: resources/v5_baits.bed
    - step: ngs_data_qc
      name: coverage_qc
      config:
        depends_on:
          alignments: mapping
          reference: genome
          kits: kits
        tool: alfred
        alfred: {}
"""

from typing import Any

from snakemake.io import directory, touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import NgsDataQc as NgsDataQcConfigModel
from .model import Tool


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

    wrappers = {"run": "fastqc"}
    result_keys = ("fastqc_done",)

    report_outputs = {
        "run": {"html": ReportOutput("FastQC", caption="fastqc.rst", htmlindex="index.html")}
    }

    default_resource_usage = ResourceUsage(threads=1, mem="4GB", runtime="4h")
    library_requirement = "extraction_type in ['dna', 'rna']"

    def _get_input_files_run(self, wildcards):
        return reads_input_files(self.parent, wildcards.library_name)

    @dictify
    def _get_output_files_run(self):
        yield "fastqc_done", touch("work/{library_name}/report/.done")
        # One report per FASTQ file, so a directory with an index of them
        yield "html", directory("work/{library_name}/report/html")

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards):
        return {"num_threads": 1} | reads_params(self.parent, wildcards.library_name)


class PicardStepPart(BaseStepPart):
    """Collect Picard metrics"""

    name = "picard"
    actions = ("prepare", "metrics")
    wrappers = {"metrics": "picard/metrics", "prepare": "picard/prepare"}

    #: ``prepare`` converts the BED files of the kit to interval lists
    extra_arguments_action = "metrics"

    @property
    def result_keys(self) -> dict[str, tuple[str, ...]]:
        """The metrics files of all programs"""
        return {"metrics": tuple(self.get_output_files("metrics"))}

    @property
    def report_outputs(self):
        programs = self.config.picard.programs
        metrics = {pgm: ReportOutput(f"Picard {pgm}", caption="picard.rst") for pgm in programs}
        return {"metrics": metrics}

    intermediate_outputs = {"prepare": ("baits", "targets")}
    resource_usage = {"metrics": ResourceUsage(threads=1, runtime="24h", mem="64GB")}
    library_requirement = "extraction_type == 'dna'"

    @dictify
    def _get_input_files_prepare(self, wildcards):
        yield from self.parent.reference_inputs().items()
        kit = self.parent.get_upstream_paths("kits", library_kit=wildcards.library_kit)
        yield "baits", kit.baits
        yield "targets", kit.targets

    @dictify
    def _get_output_files_prepare(self):
        yield "baits", "work/{library_kit}/out/{library_kit}.baits.interval_list"
        yield "targets", "work/{library_kit}/out/{library_kit}.targets.interval_list"

    def _get_log_file_prepare(self):
        return {"log": "work/{library_kit}/log/{library_kit}.prepare.log"}

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


class NgsChewStepPart(BaseStepPart):
    """Fingerprints with ngs-chew, compared across the libraries of the task"""

    name = "ngs_chew"
    actions = ("fingerprint", "compare")
    wrappers = {"compare": "ngs_chew/compare", "fingerprint": "ngs_chew/fingerprint"}

    #: ``compare`` runs ngs-chew's comparison, statistics and plots, which take no
    #: ``extra_arguments``
    extra_arguments_action = "fingerprint"

    result_keys = {
        "fingerprint": ("npz", "npz_md5"),
        "compare": tuple(table + suffix for table in COMPARE_TABLES for suffix in ("", "_md5")),
    }

    report_outputs = {
        "compare": {
            "relatedness_html": ReportOutput("ngs-chew", caption="ngs_chew_relatedness.rst"),
            "var_het_html": ReportOutput("ngs-chew", caption="ngs_chew_var_het.rst"),
        }
    }

    resource_usage = {"fingerprint": ResourceUsage(threads=1, runtime="4h", mem="2GB")}
    library_requirement = "extraction_type == 'dna'"

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

    def _get_log_file_fingerprint(self):
        return {"log": "work/{library_name}/log/{library_name}.fingerprint.log"}

    def _get_params_fingerprint(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}  # get_params adds the extra_arguments

    def _get_input_files_compare(self, wildcards):
        return {"npz": [f"work/{lib}/out/{lib}.npz" for lib in self.parent.output_entities]}

    @dictify
    def _get_output_files_compare(self):
        for table in COMPARE_TABLES:
            yield table, f"work/cohort/out/cohort.{table}.txt"
            yield f"{table}_md5", f"work/cohort/out/cohort.{table}.txt.md5"
        yield "relatedness_html", "work/cohort/report/cohort.relatedness.html"
        yield "var_het_html", "work/cohort/report/cohort.var_het.html"

    def _get_log_file_compare(self):
        return {"log": "work/cohort/log/cohort.compare.log"}


class AlfredStepPart(BaseStepPart):
    """Coverage of the targets of each library's kit with ``alfred qc``"""

    name = "alfred"
    actions = ("run",)
    wrappers = {"run": "alfred/qc"}
    default_resource_usage = ResourceUsage(threads=2, runtime="4h", mem="20GB")
    library_requirement = "extraction_type == 'dna'"

    def get_result_files(self) -> list[str]:
        """The report of each library with a kit"""
        json = "work/{library_name}/out/{library_name}.target_coverage.json.gz"
        libraries = self.parent.kit_libraries()
        return [
            path.format(library_name=lib) for lib in libraries for path in (json, json + ".md5")
        ]

    @dictify
    def _get_input_files_run(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield from self.parent.reference_inputs().items()
        yield "targets_bed", self.parent.library_kit_paths("kits", wildcards.library_name).targets

    def _get_output_files_run(self):
        json = "work/{library_name}/out/{library_name}.target_coverage.json.gz"
        return {"json": json, "json_md5": json + ".md5"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"ignore_read_groups": self.config.alfred.ignore_read_groups}


class PicardSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the Picard metrics of all DNA libraries"""

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "Picard metrics", caption="picard_summary.rst", htmlindex="index.html"
            )
        }
    }

    source = ("picard", "metrics")

    @property
    def source_keys(self) -> dict[str, str]:
        """One kind of file per program, the metrics file that the program writes first"""
        return {program: program for program in self.config.picard.programs}

    format = "picard"


class NgsChewSummaryStepPart(QcSummaryStepPart):
    """The fingerprint checks: the sheet's sex and pedigree against the ngs-chew tables"""

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "Fingerprint checks", caption="ngs_chew_summary.rst", htmlindex="index.html"
            )
        }
    }

    source = ("ngs_chew", "compare")
    source_keys = {table: table for table in COMPARE_TABLES}
    format = "ngs_chew"

    def entities(self) -> list[str]:
        """The one comparison of the cohort."""
        return ["cohort"]

    def summary_options(self) -> dict[str, Any]:
        """The sheet facts the checks compare with:
        each DNA library's donor, pedigree and sex,
        and every donor of the pedigrees with its parents and library (``""`` if unsequenced)."""
        libraries = []
        for library in self.parent.output_entities:
            row = self.parent.library_row(library)
            libraries.append(
                {
                    "library": library,
                    "donor": row["donor_name"],
                    "pedigree": row["related_group"],
                    "sex": row["sex"],
                }
            )
        donors = []
        for pedigree in dict.fromkeys(entry["pedigree"] for entry in libraries):
            for row in self.parent.pedigree_donors(pedigree, include_background=True).itertuples(
                index=False
            ):
                donors.append(
                    {
                        "donor": row.donor_name,
                        "pedigree": pedigree,
                        "father": row.father_name,
                        "mother": row.mother_name,
                        "sex": row.sex,
                        "library": row.library_name,
                    }
                )
        return {"libraries": libraries, "donors": donors}


class AlfredSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the ``alfred qc`` reports of all libraries with a kit"""

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "Target coverage", caption="alfred_summary.rst", htmlindex="index.html"
            )
        }
    }

    source = ("alfred", "run")
    source_keys = {"json": "json"}
    format = "alfred"

    def entities(self) -> list[str]:
        return self.parent.kit_libraries()


class NgsDataQcWorkflow(BaseStep):
    """Perform NGS raw data QC"""

    name = "ngs_data_qc"

    #: QC data tagged with the tool, so that a consumer can require one
    produces = [DataSignature(DataType.QC, frozenset({tool})) for tool in Tool]

    config_model_class = NgsDataQcConfigModel

    tool_parts = {
        Tool.fastqc: FastQcReportStepPart,
        Tool.picard: PicardStepPart,
        Tool.ngs_chew: NgsChewStepPart,
        Tool.alfred: AlfredStepPart,
    }

    #: Every library, not only the tumors of cancer sheets
    library_default = ""

    #: The QC summary of the files of each tool, registered after the tool's part
    summary_parts = {
        Tool.picard: PicardSummaryStepPart,
        Tool.ngs_chew: NgsChewSummaryStepPart,
        Tool.alfred: AlfredSummaryStepPart,
    }

    @classmethod
    def task_library_requirement(cls, config):
        """The libraries that the tool runs on."""
        return cls.tool_parts[config.tool].library_requirement

    @classmethod
    def task_produces(cls, config, upstream):
        """QC data tagged with the tool (``varfish_export`` reads alfred's)."""
        return (DataSignature(DataType.QC, frozenset({str(config.tool)})),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local NGS QC output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        if config.tool == "ngs_chew":
            return {"npz": f"work/{lib}/out/{lib}.npz"}
        if config.tool == "alfred":
            return {"alfred_qc": f"work/{lib}/out/{lib}.target_coverage.json.gz"}
        return {"done": f"work/{lib}/report/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        if summary := self.summary_parts.get(self.config.tool):
            self.register_sub_step_classes((summary,))
        if self.config.tool == "picard" and "CollectHsMetrics" in self.config.picard.programs:
            self.check_library_kits("kits", self.output_entities)
        if self.config.tool == "alfred":
            self.check_library_kits("kits", self.output_entities, require=False)

    def kit_libraries(self) -> list[str]:
        """Return the DNA libraries with a kit, which ``alfred`` runs on."""
        return [name for name in self.output_entities if name in self.library_kits]
