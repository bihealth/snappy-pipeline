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
    mantis_msi2 = "mantis_msi2"


class SomaticMsiCallingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    microsatellites: Annotated[
        str,
        DataSignature(DataType.REGIONS, frozenset({"microsatellites"})),
        ExpectedPathSchema(ExpectedRegions),
    ]
    """Task with the BED file of the microsatellite loci, such as ``appData/hg38/GRCh38.d1.vd1.all_loci.bed`` of MANTIS"""


class SomaticMsiCalling(SnappyStepModel):
    depends_on: SomaticMsiCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]
