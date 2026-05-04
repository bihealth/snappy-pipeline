import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


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


class HomologousRecombinationDeficiencyDependsOn(SnappyModel):
    cnv_calling: str = "cnv_calling"


class HomologousRecombinationDeficiency(SnappyStepModel):
    depends_on: HomologousRecombinationDeficiencyDependsOn = Field(
        default_factory=HomologousRecombinationDeficiencyDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.scarHRD)]

    scarHRD: ScarHRD | None = None
