# -*- coding: utf-8 -*-
"""Gene expression quantification and RNA QC

The step quantifies gene expression of RNA libraries or computes RNA QC metrics.
Each task runs one tool, selected by the ``tool`` field.
The step processes RNA libraries only.
A task takes all RNA libraries unless it sets ``library_selection``.

Quantification:

- ``featurecounts``: counts reads per gene from the BAM files, with the GTF of ``depends_on.features``.
  Its result file is a table of counts and gene lengths.
- ``salmon``: estimates expression from the FASTQ files, without prior mapping.
  It needs ``depends_on.reads``, ``depends_on.features`` and a salmon index in ``depends_on.index``.
  It needs a prebuilt index and a table that maps transcripts to genes from ``depends_on.index``.

QC and strandedness, all on the BAM files:

- ``strandedness``: infers the strandedness of each library with RSeQC ``infer_experiment.py``.
  It needs a 6-column BED file of transcripts and their strand from ``depends_on.exons``.
- ``dupradar``, ``duplication``, ``rnaseqc`` and ``stats``: estimate library complexity and duplication rate, and compute coverage and alignment metrics.
  ``dupradar`` needs BAM files with marked duplicates and reads the GTF of ``depends_on.features``.
- ``rnaseqc`` reads a GTF prepared for RNA-SeQC 1.1.8 from ``depends_on.features``, a task tagged ``rnaseqc``.
  Every record of this GTF has a ``transcript_id``, so it has no ``gene`` records.
  The GENCODE attribute ``transcript_type`` marks rRNA.

``featurecounts`` and ``dupradar`` need the strandedness of each library.
The QC tools ``duplication``, ``rnaseqc`` and ``stats`` have no strand option.
The strandedness comes from a task of this step with ``tool: strandedness``, named in ``depends_on.strandedness``.
The field ``strand`` forces the strandedness of featurecounts and dupradar instead: 0 unstranded, 1 forward, 2 reverse.
A task with ``tool: strandedness`` forces its decision with ``strandedness.strand``, as ``ngs_mapping`` does with STAR.

The example maps RNA reads with STAR, infers the strandedness, and counts with featurecounts:

.. code-block:: yaml

    tasks:
      - name: genome
        step: external_data
        config:
          produces: {type: raw, tags: [reference, dna]}
          files: {fasta: resources/genome.fa}
      - name: annotation
        step: external_data
        config:
          produces: {type: raw, tags: [features]}
          files: {gtf: resources/genes.gtf}
      - name: exons
        step: external_data
        config:
          produces: {type: regions, tags: [exons]}
          files: {bed: resources/exons.bed}
      - name: star_index
        step: external_data
        config:
          produces: {type: index, tags: [star, rna]}
          files: {index: resources/star_index}
      - name: reads
        step: external_data
        config:
          produces: {type: raw, tags: [fastq]}
          library_files: [reads.tsv]
      - name: mapping
        step: ngs_mapping
        config:
          depends_on:
            reads: reads
            reference: genome
            features: annotation
            index: star_index
            exons: exons
          tool: star
          star: {}
      - name: strandedness
        step: gene_expression_quantification
        config:
          depends_on:
            alignments: mapping
            exons: exons
          tool: strandedness
          strandedness: {}
      - name: counts
        step: gene_expression_quantification
        config:
          depends_on:
            alignments: mapping
            strandedness: strandedness
            features: annotation
          tool: featurecounts
          featurecounts: {}

Outputs
-------

For each RNA library ``{library_name}``, the task links the following result files below ``tasks/<task>/output/``.
Each file has an ``.md5`` file.
The files are in ``{library_name}/out/`` and start with ``{library_name}``:

- ``featurecounts``: ``.tsv`` and ``.tsv.summary``
- ``salmon``: ``.gene.sf`` and ``.transcript.sf``, ``.lib_format_counts.json`` (the library type that Salmon detected) and the directory ``.aux`` (its auxiliary files, such as the mapping rate in ``meta_info.json``)
- ``strandedness``: ``.tsv`` and ``.decision``
- ``dupradar``: ``.tsv``
- ``duplication``: ``.seq.DupRate.xls`` and ``.pos.DupRate.xls``
- ``rnaseqc``: ``.metrics.tsv``, and the files ``.meanCoverage_<low|medium|high>.txt``, ``.meanCoverageNorm_<low|medium|high>.txt`` and ``.gapLengthHist_<low|medium|high>.txt``
- ``stats``: ``.read_alignment_report.tsv``

Downstream tasks read the key ``tsv`` of ``featurecounts`` (data type ``expression``), the keys ``gene_sf`` and ``transcript_sf`` of ``salmon``, and the key ``decision`` of ``strandedness``.

Notes
-----

``depends_on.index`` names an ``external_data`` task with ``{type: index, tags: [salmon, rna]}`` and the keys ``index`` and ``transcript_to_gene``.
The key ``index`` is the directory of the prebuilt index, with files such as ``sa.bin``, ``rsd.bin`` and ``txpInfo.bin``.
"""

from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params
from snappy_pipeline.workflows.common.strandedness import strandedness_params
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import ExpectedStrandedness, Tool
from .model import GeneExpressionQuantification as GeneExpressionQuantificationConfigModel

#: Signature of the strandedness decision of a ``tool: strandedness`` task
STRANDEDNESS_SIGNATURE = DataSignature(DataType.QC, frozenset({"strandedness"}))
#: Signature of the expression values of ``featurecounts`` and ``salmon`` tasks
EXPRESSION_SIGNATURE = DataSignature(DataType.EXPRESSION, frozenset({"rna"}))
#: Signature of the QC tools (dupradar, duplication, rnaseqc, stats)
RNA_QC_SIGNATURE = DataSignature(DataType.QC, frozenset({"rna"}))
#: The GTF that ``rnaseqc`` reads from ``depends_on.features``, prepared for RNA-SeQC
RNASEQC_FEATURES_SIGNATURE = DataSignature(DataType.RAW, frozenset({"features", "rnaseqc"}))

#: The output files of each library, below ``work/{library_name}/out/``
OUT_PATH = "work/{library_name}/out/{library_name}"

# Extensions
EXTENSIONS = {
    "featurecounts": {
        "tsv": ".tsv",
        "tsv_md5": ".tsv.md5",
        "summary": ".tsv.summary",
        "summary_md5": ".tsv.summary.md5",
    },
    "strandedness": {
        "tsv": ".tsv",
        "tsv_md5": ".tsv.md5",
        "decision": ".decision",
        "decision_md5": ".decision.md5",
    },
    "duplication": {
        "dupl_seq": ".seq.DupRate.xls",
        "dupl_seq_md5": ".seq.DupRate.xls.md5",
        "dupl_pos": ".pos.DupRate.xls",
        "dupl_pos_md5": ".pos.DupRate.xls.md5",
    },
    "dupradar": {"dupradar": ".tsv", "dupradar_md5": ".tsv.md5"},
    "rnaseqc": {
        "rnaseqc_metrics": ".metrics.tsv",
        "rnaseqc_metrics_md5": ".metrics.tsv.md5",
        "rnaseqc_meancov_low": ".meanCoverage_low.txt",
        "rnaseqc_meancov_low_md5": ".meanCoverage_low.txt.md5",
        "rnaseqc_meancov_medium": ".meanCoverage_medium.txt",
        "rnaseqc_meancov_medium_md5": ".meanCoverage_medium.txt.md5",
        "rnaseqc_meancov_high": ".meanCoverage_high.txt",
        "rnaseqc_meancov_high_md5": ".meanCoverage_high.txt.md5",
        "rnaseqc_meannorm_low": ".meanCoverageNorm_low.txt",
        "rnaseqc_meannorm_low_md5": ".meanCoverageNorm_low.txt.md5",
        "rnaseqc_meannorm_medium": ".meanCoverageNorm_medium.txt",
        "rnaseqc_meannorm_medium_md5": ".meanCoverageNorm_medium.txt.md5",
        "rnaseqc_meannorm_high": ".meanCoverageNorm_high.txt",
        "rnaseqc_meannorm_high_md5": ".meanCoverageNorm_high.txt.md5",
        "rnaseqc_gaplen_low": ".gapLengthHist_low.txt",
        "rnaseqc_gaplen_low_md5": ".gapLengthHist_low.txt.md5",
        "rnaseqc_gaplen_medium": ".gapLengthHist_medium.txt",
        "rnaseqc_gaplen_medium_md5": ".gapLengthHist_medium.txt.md5",
        "rnaseqc_gaplen_high": ".gapLengthHist_high.txt",
        "rnaseqc_gaplen_high_md5": ".gapLengthHist_high.txt.md5",
    },
    "stats": {"stats": ".read_alignment_report.tsv", "stats_md5": ".read_alignment_report.tsv.md5"},
    "salmon": {
        "gene_sf": ".gene.sf",
        "gene_sf_md5": ".gene.sf.md5",
        "transcript_sf": ".transcript.sf",
        "transcript_sf_md5": ".transcript.sf.md5",
        "lib_format_counts": ".lib_format_counts.json",
        "lib_format_counts_md5": ".lib_format_counts.json.md5",
        # Salmon's auxiliary directory: mapping rate (meta_info.json), fragment lengths, bias
        "aux": ".aux",
    },
}


class GeneExpressionQuantificationStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Class available actions
    actions = ("run",)

    @property
    def result_keys(self) -> tuple[str, ...]:
        """All outputs of the tool"""
        return tuple(EXTENSIONS[self.name])

    #: Whether the tool reads the strandedness decision of the ``strandedness`` dependency.
    needs_strandedness = False

    def _get_input_files_run(self, wildcards: Wildcards):
        """Resolve alignment inputs through the typed upstream contract broker."""
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        input_files = {
            "bam": alignments.bam,
            "bai": alignments.bai,
        }
        if self.needs_strandedness:
            strandedness: ExpectedStrandedness = self.parent.get_upstream_paths(
                "strandedness", library_name=wildcards.library_name
            )
            input_files["decision"] = strandedness.decision
        return input_files

    def _get_output_files_run(self):
        return {key: OUT_PATH + ext for key, ext in EXTENSIONS[self.name].items()}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"strand": self.config.strand}


class StrandednessStepPart(GeneExpressionQuantificationStepPart):
    """The strand of each library with RSeQC ``infer_experiment``"""

    #: Step name
    name = "strandedness"

    wrappers = {"run": "rseqc"}
    default_resource_usage = ResourceUsage(threads=2, runtime="12h", mem="6700MB")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "exon_bed", self.parent.get_upstream_paths("exons").bed

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return strandedness_params(self.config.strandedness, wildcards.library_name)


class FeatureCountsStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "featurecounts"

    wrappers = {"run": "featurecounts"}
    default_resource_usage = ResourceUsage(threads=2, runtime="1d", mem="6700MB")
    needs_strandedness = True

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "features", self.parent.get_upstream_paths("features").gtf

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        section = self.config.featurecounts
        return super()._get_params_run(wildcards) | {
            "feature_type": section.feature_type,
            "gene_attribute": section.gene_attribute,
            "exclude_flags": section.exclude_flags,
            "min_mapq": section.min_mapq,
        }


class QCStepPartDuplication(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "duplication"
    wrappers = {"run": "rnaqc/duplication"}
    default_resource_usage = ResourceUsage(threads=1, runtime="3d", mem="127GB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}


class QCStepPartDupradar(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "dupradar"
    wrappers = {"run": "rnaqc/dupradar"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.dupradar.num_threads, runtime="4d", mem="6700MB")

    needs_strandedness = True

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "features", self.parent.get_upstream_paths("features").gtf

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return super()._get_params_run(wildcards) | {
            "num_threads": self.config.dupradar.num_threads,
        }


class QCStepPartRnaseqc(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "rnaseqc"
    wrappers = {"run": "rnaqc/rnaseqc"}
    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="16GB")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        features = self.parent.get_upstream_paths("features", signature=RNASEQC_FEATURES_SIGNATURE)
        yield "features", features.gtf

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}


class QCStepPartStats(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "stats"
    wrappers = {"run": "rnaqc/stats"}
    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="4GB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}


class SalmonStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification for raw data using salmon"""

    #: Step name
    name = "salmon"

    wrappers = {"run": "salmon"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.salmon.num_threads, runtime="4h", mem="32GB")

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        yield from reads_input_files(self.parent, wildcards.library_name).items()
        yield "features", self.parent.get_upstream_paths("features").gtf
        index = self.parent.get_upstream_paths("index")
        yield "index", index.index
        yield "transcript_to_gene", index.transcript_to_gene

    def _get_output_files_run(self):
        outputs = super()._get_output_files_run()
        outputs["aux"] = directory(outputs["aux"])
        return outputs

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result |= self.config.salmon.model_dump(by_alias=True)
        result["strand"] = self.config.strand
        return result


class GeneExpressionQuantificationWorkflow(BaseStep):
    """Perform gene expression quantification"""

    #: Workflow name
    name = "gene_expression_quantification"

    produces = [EXPRESSION_SIGNATURE, STRANDEDNESS_SIGNATURE, RNA_QC_SIGNATURE]
    config_model_class = GeneExpressionQuantificationConfigModel

    tool_parts = {
        Tool.strandedness: StrandednessStepPart,
        Tool.featurecounts: FeatureCountsStepPart,
        Tool.duplication: QCStepPartDuplication,
        Tool.dupradar: QCStepPartDupradar,
        Tool.rnaseqc: QCStepPartRnaseqc,
        Tool.stats: QCStepPartStats,
        Tool.salmon: SalmonStepPart,
    }

    library_requirement = "extraction_type == 'rna'"
    library_default = ""

    @classmethod
    def task_produces(cls, config, upstream):
        """Expression for featurecounts and salmon, the strandedness decision, or RNA QC."""
        match config.tool:
            case Tool.featurecounts | Tool.salmon:
                return (EXPRESSION_SIGNATURE,)
            case Tool.strandedness:
                return (STRANDEDNESS_SIGNATURE,)
        return (RNA_QC_SIGNATURE,)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return the local paths of the outputs of the task's tool, without the md5 files."""
        prefix = OUT_PATH.format(library_name=kwargs.get("library_name", "{library_name}"))
        return {
            key: prefix + ext
            for key, ext in EXTENSIONS[config.tool].items()
            if not key.endswith("_md5")
        }
