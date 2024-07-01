import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    manta = "manta"
    delly2 = "delly2"


class Manta(SnappyModel):
    pass


class Delly2(SnappyModel):
    path_exclude_tsv: str | None = None
    max_threads: int = 16


class SomaticWgsSvCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.manta], min_length=1)]

    manta: Manta | None = None

    delly2: Delly2 | None = None
