import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    fastqc = "fastqc"
    picard = "picard"


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


class Picard(SnappyModel):
    path_to_baits: str = ""
    """Required when CollectHsMetrics is among the programs"""

    path_to_targets: str = ""
    """When missing, same as baits"""

    bait_name: str = ""
    """Exon enrichment kit name (optional)"""

    programs: Annotated[list[PicardProgram], Field(min_length=1)]

    @model_validator(mode="after")
    def ensure_baits_when_required(self):
        if PicardProgram.CollectHsMetrics in self.programs and not self.path_to_baits:
            raise ValueError(
                "Path to baits is required when CollectHsMetrics is among the programs"
            )
        return self


class Fastqc(SnappyModel):
    pass


class NgsDataQcDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    link_in: str | None = None
    """Optional: name of the ``link_in`` task to use as the preprocessed FASTQ source."""


class NgsDataQc(SnappyStepModel):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: NgsDataQcDependsOn = Field(default_factory=NgsDataQcDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.fastqc)]

    picard: Picard | None = None

    fastqc: Fastqc | None = None  # TODO fastqc has no configuration options in the DEFAULT_CONFIG?
