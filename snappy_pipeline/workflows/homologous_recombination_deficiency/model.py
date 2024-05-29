import enum
from typing import Annotated

from pydantic import model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    scarHRD = "scarHRD"


class GenomeName(enum.StrEnum):
    grch37 = "grch37"
    grch38 = "grch38"
    mouse = "mouse"


class ScarHRD(SnappyModel):
    genome_name: GenomeName = GenomeName.grch37

    chr_prefix: bool = False

    length: int = 50
    """Wiggle track for GC reference file"""


class HomologousRecombinationDeficiency(SnappyStepModel):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.scarHRD], min_length=1)]

    path_cnv_calling: str

    scarHRD: ScarHRD | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
