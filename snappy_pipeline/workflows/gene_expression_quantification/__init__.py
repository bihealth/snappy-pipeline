# -*- coding: utf-8 -*-
"""Implementation of the ``gene_expression_quantification`` step

This step takes as input alignments from ``ngs_mapping`` and performs gene
expression quantification. Further, the tool dupradar is run to help with
estimating library complexity, i.e. the number of PCR duplicates.
(This needs the bam files to be marked for duplicates, e.g. with samblaster.)

The output are tables of gene counts and gene lengths for each library for
a given annotation.

==========
Step Input
==========

The gene expression quantification step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step. It will use RNA-seq experiments only.

Additionally, salmon can be used to estimate expression directly from the FASTQs with no need
for prior mapping.

The tools ``featurecounts``, ``dupradar``, ``duplication``, ``rnaseqc`` and ``stats`` also need
the strandedness of each library. It comes from a separate task of this step with
``tool: strandedness``, named in ``depends_on.strandedness``:

.. code-block:: yaml

    tasks:
      - step: gene_expression_quantification
        name: strandedness
        config:
          depends_on:
            alignments: mapping
          tool: strandedness
          strandedness:
            path_exon_bed: /path/to/exons.bed
      - step: gene_expression_quantification
        name: counts
        config:
          depends_on:
            alignments: mapping
            strandedness: strandedness
          tool: featurecounts

===========
Step Output
===========

For each RNA NGS library with name ``lib_name``/key ``lib_pk``, the
pipeline step will create a directory ``output/{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting TSV files.

- ``{lib_name}-{lib_pk}.tsv``

=====================
Default Configuration
=====================

Featurecounts needs a GTF with the gene model that will be used to count reads. By default,
strandedness is set to 0, i.e. unstranded. If you want stranded output, see the featurecounts
manual and set it appropriately (1 or 2, depending on the protocoll).

RSeQC's infer_experiment.py can be used to infer the strandedness. (You might want to run it
separately before featurecounts). This requires a 6-column bed file, which holds the positions
of transcripts and the strand they are originating from.

Salmon needs an pre-build index at the moment, the path should point to the directory, which
contains the needed files (e.g. sa.bin, rsd.bin, txpInfo.bin, etc).
Additionally, one can provide a gtf for the mapping between transcripts and genes.
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
from .model import Salmon as SalmonConfigModel

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

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: SalmonConfigModel = self.config.salmon

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        if self.config.tool != self.name:
            return
        yield from reads_input_files(self.parent, wildcards.library_name).items()
        yield "features", self.parent.get_upstream_paths("features").gtf
        if self.cfg and self.cfg.path_index:
            yield "index", self.cfg.path_index
        if self.cfg and self.cfg.path_transcript_to_gene:
            yield "transcript_to_gene", self.cfg.path_transcript_to_gene

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result |= self.config.salmon.model_dump(
            by_alias=True, exclude={"path_index", "path_transcript_to_gene"}
        )
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
        yield "exon_bed", self.config.strandedness.path_exon_bed

    default_resource_usage = ResourceUsage(threads=2, runtime="12h", mem="6700MB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if self.config.tool != self.name:
            return super()._get_params_run(wildcards)
        config = self.config.strandedness.model_dump(by_alias=True, exclude={"path_exon_bed"}) | {
            "strand": self.config.strand
        }
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
        yield "reference", self.parent.get_upstream_paths("reference").fasta
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

        df = self.build_library_dataframe()
        if df.empty:
            return []

        all_fns = []
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type.lower() == "rna":
                all_fns += expand(
                    "output/{library_name}/out/{library_name}{ext}",
                    library_name=[library_name],
                    ext=exts,
                )

        return all_fns
