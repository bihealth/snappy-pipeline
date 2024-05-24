import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel, EnumField


class Tool(enum.StrEnum):
    mantis = "mantis"


class SomaticMsiCalling(SnappyStepModel):
    path_ngs_mapping: str

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mantis])]

    loci_bed: Annotated[
        str, Field(examples=["/fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed"])
    ]
