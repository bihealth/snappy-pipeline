import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Exons,
    ExpectedPathSchema,
    Features,
    Reference,
)
from snappy_pipeline.workflows.common.strandedness import Strand, Strandedness
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Featurecounts(ToolModel):
    feature_type: str = "exon"
    """Feature type (GTF column 3) to count (``featureCounts -t``)"""

    gene_attribute: str = "gene_id"
    """GTF attribute that groups the features into genes (``featureCounts -g``)"""

    exclude_flags: int = 260
    """SAM flags of the reads to leave out: unmapped and secondary alignments (``samtools view -F``)"""

    min_mapq: int = 255
    """Minimum mapping quality of the reads (``samtools view -q``); 255 keeps the uniquely mapped reads of STAR"""


class RnaSeqC(ToolModel):
    pass


class DupRadar(SnappyModel):
    num_threads: int = 8


class Salmon(ToolModel):
    salmon_params: str = " --gcBias --validateMappings"
    num_threads: int = 16


class Duplication(ToolModel):
    pass


class Stats(SnappyModel):
    pass


class ExpectedExpression(SnappyModel):
    """Consumer-driven contract: expected output keys from gene_expression_quantification."""

    tsv: str


class ExpectedSalmonIndex(SnappyModel):
    """Consumer-driven contract: the salmon index and its transcript to gene table."""

    index: str
    """Directory of the prebuilt index"""

    transcript_to_gene: str


class Tool(enum.StrEnum):
    strandedness = "strandedness"
    featurecounts = "featurecounts"
    dupradar = "dupradar"
    duplication = "duplication"
    rnaseqc = "rnaseqc"
    salmon = "salmon"
    stats = "stats"


class ExpectedStrandedness(SnappyModel):
    """Consumer-driven contract: the strandedness decision (JSON) of an RNA library."""

    decision: str


class GeneExpressionQuantificationDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedAlignments),
        RequiredBy(
            Tool.strandedness,
            Tool.featurecounts,
            Tool.dupradar,
            Tool.duplication,
            Tool.rnaseqc,
            Tool.stats,
        ),
    ] = ""
    """Alignments of RNA libraries, from an ``ngs_mapping`` task"""

    reads: Annotated[str, DataSignature(DataType.RAW), RequiredBy(Tool.salmon)] = ""
    """FASTQ source: ``data_sets``, an ``external_data`` task or an ``adapter_trimming`` task"""

    strandedness: Annotated[
        str,
        DataSignature(DataType.QC, frozenset({"strandedness"})),
        ExpectedPathSchema(ExpectedStrandedness),
        RequiredBy(Tool.featurecounts, Tool.dupradar),
    ] = ""
    """Task of this step with ``tool: strandedness``"""

    reference: Annotated[Reference, RequiredBy(Tool.rnaseqc)] = ""

    features: Annotated[
        Features, RequiredBy(Tool.featurecounts, Tool.dupradar, Tool.rnaseqc, Tool.salmon)
    ] = ""
    """The gene annotation; for ``tool: rnaseqc``, a GTF that is also tagged ``rnaseqc``: RNA-SeQC
    1.1.8 needs a ``transcript_id`` in every record (no ``gene`` records) and finds rRNA by the
    GENCODE attribute ``transcript_type``"""

    exons: Annotated[Exons, RequiredBy(Tool.strandedness)] = ""
    """Task with a BED file of transcripts with their strand in column 6"""

    index: Annotated[
        str,
        DataSignature(DataType.INDEX, frozenset({"salmon", "rna"})),
        ExpectedPathSchema(ExpectedSalmonIndex),
        RequiredBy(Tool.salmon),
    ] = ""
    """An ``external_data`` task with the prebuilt index and the table that maps transcripts to
    genes"""


class GeneExpressionQuantification(SnappyStepModel):
    depends_on: GeneExpressionQuantificationDependsOn = Field(
        default_factory=GeneExpressionQuantificationDependsOn
    )

    tool: Tool

    strand: Strand = Strand.infer
    """Use 0, 1 or 2 to force unstranded, forward or reverse strand for ``featurecounts``, ``dupradar`` and ``salmon``.
    Use -1 to take the decision of ``depends_on.strandedness`` (``salmon``: to guess)."""

    featurecounts: Featurecounts | None = None

    strandedness: Strandedness | None = None

    rnaseqc: RnaSeqC | None = None

    dupradar: DupRadar | None = None

    duplication: Duplication | None = None

    stats: Stats | None = None

    salmon: Salmon | None = None

    @model_validator(mode="after")
    def _strandedness_forces_its_own_strand(self):
        """``tool: strandedness`` forces the strand with ``strandedness.strand``, like ``ngs_mapping``"""
        if self.tool == Tool.strandedness and self.strand != Strand.infer:
            raise ValueError(
                "tool=strandedness takes a forced strand from strandedness.strand, not strand"
            )
        return self
