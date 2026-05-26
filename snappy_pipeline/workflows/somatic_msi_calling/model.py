import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    mantis_msi2 = "mantis_msi2"


class SomaticMsiCallingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class SomaticMsiCalling(SnappyStepModel):
    depends_on: SomaticMsiCallingDependsOn = Field(default_factory=SomaticMsiCallingDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.mantis_msi2)]

    loci_bed: Annotated[
        str,
        Field(
            examples=[
                "/fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed",
                "/fast/work/groups/cubi/projects/biotools/Mantis/appData/hg38/GRCh38.d1.vd1.all_loci.bed",
            ]
        ),
    ]
