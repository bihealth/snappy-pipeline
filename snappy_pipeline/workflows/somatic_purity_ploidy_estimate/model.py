import enum
from typing import Annotated


from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedRegions,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    ascat = "ascat"


class SomaticPurityPloidyEstimateDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    snp_loci: Annotated[
        str,
        DataSignature(DataType.REGIONS, frozenset({"snp_loci"})),
        ExpectedPathSchema(ExpectedRegions),
    ]
    """Task with the BED file of the loci for B allele frequency"""


class SomaticPurityPloidyEstimate(SnappyStepModel):
    depends_on: SomaticPurityPloidyEstimateDependsOn

    tool: Annotated[Tool, EnumField(Tool)]
