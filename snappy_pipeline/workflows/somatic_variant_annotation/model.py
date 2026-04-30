import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.annotation import Mehari, Vep


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class SomaticVariantAnnotationDependsOn(SnappyModel):
    somatic_variant: str = "somatic_variant"


class SomaticVariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    depends_on: SomaticVariantAnnotationDependsOn = Field(
        default_factory=SomaticVariantAnnotationDependsOn
    )

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    vep: Vep | None = None

    mehari: Mehari | None = None
