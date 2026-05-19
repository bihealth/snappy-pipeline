import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    optitype = "optitype"
    arcashla = "arcashla"


class Optitype(SnappyModel):
    max_reads: int = 5000
    """5000 is a suggestion by OptiType author"""

    num_mapping_threads: int = 4


class ArcasHla(SnappyModel):
    mapper: str = "star"


class HlaTypingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    link_in: str | None = None
    """Optional: name of the ``link_in`` task to use as the preprocessed FASTQ source."""


class HlaTyping(SnappyStepModel):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: HlaTypingDependsOn = Field(default_factory=HlaTypingDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.optitype)]

    optitype: Optitype = Optitype()

    arcashla: ArcasHla = ArcasHla()
