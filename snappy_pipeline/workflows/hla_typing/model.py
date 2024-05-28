import enum
from typing import Annotated, Self

from pydantic import model_validator

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


class HlaTyping(SnappyStepModel):
    path_ngs_mapping: str

    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.optitype], min_length=1)]

    optitype: Optitype | None = None

    arcashla: ArcasHla | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
