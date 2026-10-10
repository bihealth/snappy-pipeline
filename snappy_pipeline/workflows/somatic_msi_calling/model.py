import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
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


class SomaticMsiCalling(SnappyStepModel):
    depends_on: SomaticMsiCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    loci_bed: Annotated[
        str,
        Field(
            examples=[
                "/path/to/Mantis/appData/hg19/loci.bed",
                "/path/to/Mantis/appData/hg38/GRCh38.d1.vd1.all_loci.bed",
            ]
        ),
    ]
