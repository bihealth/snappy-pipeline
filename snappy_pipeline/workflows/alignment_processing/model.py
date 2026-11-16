import enum
from typing import Annotated


from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    CommonVariants,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    gatk4_bqsr = "gatk4_bqsr"


class AlignmentProcessingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    """The BAM files to process: an ``ngs_mapping`` task or another ``alignment_processing`` task"""

    reference: Reference

    common_variants: CommonVariants
    """Known variant sites, which BQSR does not count as sequencing errors (``gatk4_bqsr``)"""


class Gatk4Bqsr(ToolModel):
    """Base quality score recalibration with GATK ``BaseRecalibrator`` and ``ApplyBQSR``;
    ``extra_arguments`` go to ``BaseRecalibrator``"""


class AlignmentProcessing(SnappyStepModel):
    bam_intermediate = True

    depends_on: AlignmentProcessingDependsOn

    tool: Tool

    gatk4_bqsr: Gatk4Bqsr | None = None
