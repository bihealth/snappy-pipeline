import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel, validators
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.adapter_trimming.model import ExpectedTrimmedRawFastq
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


class GeneExpressionQuantificationDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"

    # Optional external/in-pipeline FASTQ source for salmon mode.
    link_in: Annotated[
        str,
        DataSignature(DataType.RAW),
        ExpectedPathSchema(ExpectedLinkedRawFastq),
    ] = ""
    adapter_trimming: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"trimmed"})),
        ExpectedPathSchema(ExpectedTrimmedRawFastq),
    ] = ""


class GeneExpressionQuantification(SnappyStepModel, validators.NgsMappingMixin):
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
