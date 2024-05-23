import enum
from typing import Annotated

from pydantic import model_validator

from models import SnappyStepModel, EnumField, SnappyModel


class Tool(enum.Enum):
    fastqc = "fastqc"
    picard = "picard"


class PicardProgram(enum.Enum):
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
    path_ngs_mapping: str

    path_to_baits: str = ""
    """Required when CollectHsMetrics is among the programs"""

    path_to_targets: str = ""
    """When missing, same as baits"""

    bait_name: str = ""
    """Exon enrichment kit name (optional)"""

    programs: list[PicardProgram] = []


class Fastqc(SnappyModel):
    pass


class NgsDataQc(SnappyStepModel):
    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.fastqc, Tool.picard])]

    picard: Picard | None = None

    fastqc: Fastqc | None = None  # TODO fastqc has no configuration options in the DEFAULT_CONFIG?

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self