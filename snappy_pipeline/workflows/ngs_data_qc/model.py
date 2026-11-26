import enum
from typing import Annotated

from pydantic import BaseModel, Field

from snappy_pipeline.models import (
    EnumField,
    ReadBy,
    RequiredBy,
    SnappyModel,
    SnappyStepModel,
    ToolModel,
)
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Dbsnp,
    DnaAlignments,
    Kits,
    Reference,
)


class Tool(enum.StrEnum):
    fastqc = "fastqc"
    picard = "picard"
    ngs_chew = "ngs_chew"
    alfred = "alfred"


class PicardProgram(enum.StrEnum):
    # Generic metrics
    CollectAlignmentSummaryMetrics = "CollectAlignmentSummaryMetrics"
    CollectBaseDistributionByCycle = "CollectBaseDistributionByCycle"
    CollectGcBiasMetrics = "CollectGcBiasMetrics"
    CollectInsertSizeMetrics = "CollectInsertSizeMetrics"
    CollectQualityYieldMetrics = "CollectQualityYieldMetrics"
    CollectSequencingArtifactMetrics = "CollectSequencingArtifactMetrics"
    MeanQualityByCycle = "MeanQualityByCycle"
    QualityScoreDistribution = "QualityScoreDistribution"

    # The above are grouped into
    CollectMultipleMetrics = "CollectMultipleMetrics"

    # Generic metrics not included in "CollectMultipleMetrics"
    CollectJumpingLibraryMetrics = "CollectJumpingLibraryMetrics"
    CollectOxoGMetrics = "CollectOxoGMetrics"
    EstimateLibraryComplexity = "EstimateLibraryComplexity"

    # WGS-specific metrics
    CollectRawWgsMetrics = "CollectRawWgsMetrics"
    CollectWgsMetrics = "CollectWgsMetrics"
    CollectWgsMetricsWithNonZeroCoverage = "CollectWgsMetricsWithNonZeroCoverage"

    # Other assay-specific metrics
    CollectHsMetrics = "CollectHsMetrics"
    """Whole Exome Sequencing"""

    CollectTargetedPcrMetrics = "CollectTargetedPcrMetrics"
    """Panel sequencing"""

    CollectRnaSeqMetrics = "CollectRnaSeqMetrics"
    """mRNA sequencing, not implemented yet"""

    CollectRbsMetrics = "CollectRbsMetrics"
    """bi-sulfite sequencing, not implemented yet"""


class Picard(ToolModel):
    """``extra_arguments`` go to every Picard program of ``programs``"""

    programs: Annotated[list[PicardProgram], Field(min_length=1)]


class Fastqc(ToolModel):
    pass


class Alfred(ToolModel):
    """Coverage of the targets of each library's kit with ``alfred qc``"""

    ignore_read_groups: bool = True
    """Compute one set of metrics per library, not per read group (``alfred qc --ignore``)"""


class ExpectedAlfredQc(BaseModel):
    """Consumer-driven contract: the ``alfred qc`` report of an ``ngs_data_qc`` task."""

    alfred_qc: str


class NgsChew(ToolModel):
    """``extra_arguments`` go to ``ngs-chew fingerprint``, for example
    ``["--genome-release", "GRCh38"]`` when ngs-chew cannot tell the release from the contigs;
    the reference must have all contigs of the release, since ngs-chew calls its sites on all of
    them"""


class NgsDataQcDependsOn(SnappyModel):
    alignments: Annotated[DnaAlignments, RequiredBy(Tool.picard, Tool.ngs_chew, Tool.alfred)] = ""
    """Alignments of DNA libraries"""

    reads: Annotated[str, DataSignature(DataType.RAW), RequiredBy(Tool.fastqc)] = ""
    """FASTQ source: an ``external_data`` task or an ``adapter_trimming`` task"""

    reference: Annotated[Reference, RequiredBy(Tool.picard, Tool.ngs_chew, Tool.alfred)] = ""

    dbsnp: Annotated[Dbsnp, ReadBy(Tool.picard)] = ""

    kits: Annotated[
        Kits,
        RequiredBy(Tool.alfred),
        RequiredBy(Tool.picard, when="picard.programs", contains=PicardProgram.CollectHsMetrics),
    ] = ""
    """Baits and targets of the library kits"""


class NgsDataQc(SnappyStepModel):
    """Quality control of reads and alignments"""

    depends_on: NgsDataQcDependsOn = Field(default_factory=NgsDataQcDependsOn)

    tool: Annotated[Tool, EnumField(Tool)]

    picard: Picard | None = None

    fastqc: Fastqc | None = None

    ngs_chew: NgsChew | None = None

    alfred: Alfred | None = None
