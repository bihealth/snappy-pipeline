import enum
from typing import Annotated

from pydantic import model_validator, Field, DirectoryPath

from snappy_pipeline.models import SnappyStepModel, EnumField, SnappyModel


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
    path_ngs_mapping: DirectoryPath | str

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


class NgsDataQc(SnappyStepModel):
    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.fastqc, Tool.picard], min_length=1)]

    picard: Picard | None = None

    fastqc: Fastqc | None = None  # TODO fastqc has no configuration options in the DEFAULT_CONFIG?

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
