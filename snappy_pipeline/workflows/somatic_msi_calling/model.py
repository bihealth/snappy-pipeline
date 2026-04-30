import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    mantis_msi2 = "mantis_msi2"


class SomaticMsiCallingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class SomaticMsiCalling(SnappyStepModel):
    depends_on: SomaticMsiCallingDependsOn = Field(default_factory=SomaticMsiCallingDependsOn)

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mantis_msi2], min_length=1)]

    loci_bed: Annotated[
        str,
        Field(
            examples=[
                "/fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed",
                "/fast/work/groups/cubi/projects/biotools/Mantis/appData/hg38/GRCh38.d1.vd1.all_loci.bed",
            ]
        ),
    ]
