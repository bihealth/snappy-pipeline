import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import (
    EnumField,
    SnappyModel,
    SnappyStepModel,
    ToolModel,
    validators,
)
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Dbsnp,
    ExpectedPathSchema,
    Kits,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


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


class Picard(ToolModel):
    """``extra_arguments`` go to every Picard program of ``programs``"""

    programs: Annotated[list[PicardProgram], Field(min_length=1)]


class Fastqc(ToolModel):
    pass


class NgsDataQcDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""
    """Alignments of DNA libraries, read by ``picard``"""

    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
    ] = ""
    """FASTQ source: ``data_sets``, an ``external_data`` task or an ``adapter_trimming`` task"""

    reference: Reference = ""

    dbsnp: Dbsnp = ""

    kits: Kits = ""
    """Baits and targets of the library kits, read by ``picard`` for ``CollectHsMetrics``"""


#: The ``depends_on`` fields each tool reads.
TOOL_DEPENDENCIES = {Tool.fastqc: ("reads",), Tool.picard: ("alignments", "reference")}


class NgsDataQc(SnappyStepModel):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: NgsDataQcDependsOn = Field(default_factory=NgsDataQcDependsOn)

    tool: Annotated[Tool, EnumField(Tool)]

    picard: Picard | None = None

    fastqc: Fastqc | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self

    @model_validator(mode="after")
    def ensure_kits_for_hs_metrics(self):
        if (
            self.picard
            and PicardProgram.CollectHsMetrics in self.picard.programs
            and not self.depends_on.kits
        ):
            raise ValueError("CollectHsMetrics needs depends_on.kits")
        return self
