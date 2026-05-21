import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    manta = "manta"
    delly2 = "delly2"


class Manta(SnappyModel):
    pass


class Delly2(SnappyModel):
    path_exclude_tsv: str | None = None
    max_threads: int = 16


class SomaticWgsSvCallingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class SomaticWgsSvCalling(SnappyStepModel):
    depends_on: SomaticWgsSvCallingDependsOn = Field(default_factory=SomaticWgsSvCallingDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.manta)]

    manta: Manta | None = None

    delly2: Delly2 | None = None
