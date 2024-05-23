import enum
from typing import Annotated

from pydantic import Field

from models import SnappyStepModel, EnumField, SnappyModel


class Tool(enum.Enum):
    manta = "manta"
    delly2 = "delly2"


class Manta(SnappyModel):
    pass


class Delly2(SnappyModel):
    path_exclude_tsv: str | None = None
    max_threads: int = 16


class SomaticWgsSvCalling(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.manta])]

    manta: Manta | None = None

    delly2: Delly2 | None = None
