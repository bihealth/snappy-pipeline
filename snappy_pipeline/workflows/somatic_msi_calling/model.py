import enum
from typing import Annotated

from pydantic import Field, DirectoryPath

from snappy_pipeline.models import SnappyStepModel, EnumField


class Tool(enum.StrEnum):
    mantis = "mantis"


class SomaticMsiCalling(SnappyStepModel):
    path_ngs_mapping: DirectoryPath | str

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mantis], min_length=1)]

    loci_bed: Annotated[
        str, Field(examples=["/fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed"])
    ]
