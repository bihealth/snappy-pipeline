# -*- coding: utf-8 -*-
"""Gene expression quantification and RNA QC

The step quantifies gene expression of RNA libraries or computes RNA QC metrics.
Each task runs one tool, selected by the ``tool`` field.
The step uses RNA libraries only; a task takes all RNA libraries unless it sets ``library_selection``.

Quantification:

- ``featurecounts``: counts reads per gene from the BAM files, with the GTF of ``depends_on.features``. It publishes a table of counts and gene lengths.
- ``salmon``: estimates expression from the FASTQ files, without prior mapping. It needs a prebuilt index and a table that maps transcripts to genes from ``depends_on.index``.

QC and strandedness, all on the BAM files:

- ``strandedness``: infers the strandedness of each library with RSeQC ``infer_experiment.py``. It needs a 6-column BED file of transcripts and their strand from ``depends_on.exons``.
- ``dupradar``, ``duplication``, ``rnaseqc`` and ``stats``: estimate library complexity and duplication rate, and compute coverage and alignment metrics. ``dupradar`` needs BAM files marked for duplicates. ``rnaseqc`` and ``dupradar`` read a GTF from their own config sections.

``featurecounts``, ``dupradar``, ``duplication``, ``rnaseqc`` and ``stats`` need the strandedness of each library.
It comes from a task of this step with ``tool: strandedness``, named in ``depends_on.strandedness``.
The field ``strand`` forces the strandedness of featurecounts (0 unstranded, 1 forward, 2 reverse) instead.

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
      - name: mapping
        step: ngs_mapping
        config:
          depends_on:
            reads: data_sets
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

The task publishes the following files below ``tasks/<task>/output/``, for each RNA library with name ``{library_name}``.
Each file has an ``.md5`` file.
The files are in ``{library_name}/out/`` and start with ``{library_name}``:

- ``featurecounts``: ``.tsv`` and ``.tsv.summary``
- ``salmon``: ``.gene.sf`` and ``.transcript.sf``
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

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import ExpectedStrandedness, Tool
from .model import GeneExpressionQuantification as GeneExpressionQuantificationConfigModel

#: Signature of the strandedness decision of a ``tool: strandedness`` task
STRANDEDNESS_SIGNATURE = DataSignature(DataType.QC, frozenset({"strandedness"}))
#: Signature of the expression values of ``featurecounts`` and ``salmon`` tasks
EXPRESSION_SIGNATURE = DataSignature(DataType.EXPRESSION, frozenset({"rna"}))
#: Signature of the QC tools (dupradar, duplication, rnaseqc, stats)
RNA_QC_SIGNATURE = DataSignature(DataType.QC, frozenset({"rna"}))

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
    },
}


class GeneExpressionQuantificationStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Whether the tool reads the strandedness decision of the ``strandedness`` dependency.
    needs_strandedness = False

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"

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
        if self.name != self.config.tool:
            return {}
        return {
            key: self.base_path_out.format(ext=ext) for key, ext in EXTENSIONS[self.name].items()
        }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"strand": self.config.strand}

    def _get_log_file_run(self):
        if self.name != self.config.tool:
            return {}
        return {"log": "work/{library_name}/log/{library_name}.log"}


class SalmonStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification for raw data using salmon"""

    #: Step name
    name = "salmon"

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        if self.config.tool != self.name:
            return
        yield from reads_input_files(self.parent, wildcards.library_name).items()
        yield "features", self.parent.get_upstream_paths("features").gtf
        index = self.parent.get_upstream_paths("index")
        yield "index", index.index
        yield "transcript_to_gene", index.transcript_to_gene

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result |= self.config.salmon.model_dump(by_alias=True)
        result["strand"] = self.config.strand
        return result

    default_resource_usage = ResourceUsage(threads=16, runtime="4h", mem="32GB")


class FeatureCountsStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "featurecounts"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "features", self.parent.get_upstream_paths("features").gtf

    default_resource_usage = ResourceUsage(threads=2, runtime="1d", mem="6700MB")


class StrandednessStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "strandedness"

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.tool != self.name:
            return
        yield "exon_bed", self.parent.get_upstream_paths("exons").bed

    default_resource_usage = ResourceUsage(threads=2, runtime="12h", mem="6700MB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if self.config.tool != self.name:
            return super()._get_params_run(wildcards)
        config = self.config.strandedness.model_dump(by_alias=True) | {"strand": self.config.strand}
        return {"config": config, "library_name": wildcards.library_name}


class QCStepPartDuplication(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "duplication"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, runtime="3d", mem="127GB")


class QCStepPartDupradar(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "dupradar"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.tool != self.name:
            return
        yield "dupradar_path_annotation_gtf", self.config.dupradar.dupradar_path_annotation_gtf

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if self.config.tool != self.name:
            return super()._get_params_run(wildcards)
        return super()._get_params_run(wildcards) | {
            "num_threads": self.config.dupradar.num_threads,
        }

    default_resource_usage = ResourceUsage(threads=8, runtime="4d", mem="6700MB")


class QCStepPartRnaseqc(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "rnaseqc"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.tool != self.name:
            return
        yield from self.parent.reference_inputs().items()
        yield "rnaseqc_path_annotation_gtf", self.config.rnaseqc.rnaseqc_path_annotation_gtf

    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="16GB")


class QCStepPartStats(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "stats"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="4GB")


class GeneExpressionQuantificationWorkflow(BaseStep):
    """Perform gene expression quantification"""

    #: Workflow name
    name = "gene_expression_quantification"

    config_model_class = GeneExpressionQuantificationConfigModel

    produces = [EXPRESSION_SIGNATURE, STRANDEDNESS_SIGNATURE, RNA_QC_SIGNATURE]

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

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
        """Return the local paths of the expression or strandedness outputs of the task's tool."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/{lib}/out/{lib}"
        match config.tool:
            case Tool.featurecounts:
                return {"tsv": f"{prefix}.tsv"}
            case Tool.salmon:
                return {"gene_sf": f"{prefix}.gene.sf", "transcript_sf": f"{prefix}.transcript.sf"}
            case Tool.strandedness:
                return {"decision": f"{prefix}.decision"}
            case Tool.dupradar:
                return {"tsv": f"{prefix}.tsv"}
        return {}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                StrandednessStepPart,
                FeatureCountsStepPart,
                QCStepPartDuplication,
                QCStepPartDupradar,
                QCStepPartRnaseqc,
                QCStepPartStats,
                SalmonStepPart,
            )
        )
        # Inputs are resolved via get_upstream_paths() in step parts.

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        exts = EXTENSIONS[self.config.tool].values()

        return expand(
            "output/{library_name}/out/{library_name}{ext}",
            library_name=self.output_entities,
            ext=exts,
        )
