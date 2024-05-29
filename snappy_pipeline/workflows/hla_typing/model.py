import enum
from typing import Annotated

from models.validators import NgsMappingMixin, ToolsMixin

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


class HlaTyping(SnappyStepModel, ToolsMixin, NgsMappingMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.optitype], min_length=1)]

    optitype: Optitype | None = None

    arcashla: ArcasHla | None = None
