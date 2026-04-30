import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


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


class HlaTyping(SnappyStepModel, validators.ToolsMixin, validators.NgsMappingMixin):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: HlaTypingDependsOn = Field(default_factory=HlaTypingDependsOn)

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.optitype], min_length=1)]

    optitype: Optitype = Optitype()

    arcashla: ArcasHla = ArcasHla()
