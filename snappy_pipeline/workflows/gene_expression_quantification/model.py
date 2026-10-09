import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, validators
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.link_in.model import ExpectedLinkedRawFastq
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Strand(enum.IntEnum):
    unstranded = 0
    forward = 1
    reverse = 2


class Featurecounts(SnappyModel):
    path_annotation_gtf: str


class Strandedness(SnappyModel):
    path_exon_bed: str
    """needs column 6 with strand info, e.g. CCDS/15/GRCh37/CCDS.bed"""

    threshold: float = 0.85


class RnaSeqC(SnappyModel):
    rnaseqc_path_annotation_gtf: str


class DupRadar(SnappyModel):
    dupradar_path_annotation_gtf: str
    num_threads: int = 8


class Salmon(SnappyModel):
    path_index: str
    path_transcript_to_gene: str
    salmon_params: str = " --gcBias --validateMappings"
    num_threads: int = 16


class Duplication(SnappyModel):
    pass


class Stats(SnappyModel):
    pass


class ExpectedExpression(SnappyModel):
    """Consumer-driven contract: expected output keys from gene_expression_quantification."""

    tsv: str


class Tool(enum.StrEnum):
    strandedness = "strandedness"
    featurecounts = "featurecounts"
    dupradar = "dupradar"
    duplication = "duplication"
    rnaseqc = "rnaseqc"
    salmon = "salmon"
    stats = "stats"


#: The ``depends_on`` fields each tool reads. ``strandedness`` is the decision of a
#: ``tool: strandedness`` task.
TOOL_DEPENDENCIES = {
    Tool.strandedness: ("alignments",),
    Tool.featurecounts: ("alignments", "strandedness"),
    Tool.dupradar: ("alignments", "strandedness"),
    Tool.duplication: ("alignments", "strandedness"),
    Tool.rnaseqc: ("alignments", "strandedness"),
    Tool.stats: ("alignments", "strandedness"),
    Tool.salmon: ("reads",),
}


class ExpectedStrandedness(SnappyModel):
    """Consumer-driven contract: the strandedness decision (JSON) of an RNA library."""

    decision: str


class GeneExpressionQuantificationDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""

    #: FASTQ source: a ``link_in`` task, a task whose ``output/`` holds FASTQs (such as
    #: ``adapter_trimming``), or ``data_sets`` to search the data sets' search paths.
    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
        ExpectedPathSchema(ExpectedLinkedRawFastq),
    ] = ""

    # Task of this step with tool: strandedness; required for TOOLS_NEEDING_STRANDEDNESS.
    strandedness: Annotated[
        str,
        DataSignature(DataType.QC, frozenset({"strandedness"})),
        ExpectedPathSchema(ExpectedStrandedness),
    ] = ""


class GeneExpressionQuantification(SnappyStepModel):
    depends_on: GeneExpressionQuantificationDependsOn = Field(
        default_factory=GeneExpressionQuantificationDependsOn
    )

    tool: Tool  # TODO: add default = [Tool.salmon]

    strand: Strand | int = -1  # TODO: what is this default value of -1?
    """Use 0, 1 or 2 to force unstranded, forward or reverse strand. Use -1 to guess."""

    featurecounts: Featurecounts | None = None

    strandedness: Strandedness | None = None

    rnaseqc: RnaSeqC | None = None

    dupradar: DupRadar | None = None

    duplication: Duplication | None = None

    stats: Stats | None = None

    salmon: Salmon | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
