import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class SomaticVariantAnnotationDependsOn(SnappyModel):
    somatic_variant: str = "somatic_variant"


class SomaticVariantAnnotation(SnappyStepModel):
    depends_on: SomaticVariantAnnotationDependsOn = Field(
        default_factory=SomaticVariantAnnotationDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.vep)]

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    vep: Vep | None = None

    mehari: Mehari | None = None
