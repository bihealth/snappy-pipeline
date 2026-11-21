import enum


from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    CommonVariants,
    DnaAlignments,
    Reference,
)


class Tool(enum.StrEnum):
    gatk4_bqsr = "gatk4_bqsr"


class AlignmentProcessingDependsOn(SnappyModel):
    alignments: DnaAlignments
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
