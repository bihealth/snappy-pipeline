import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyStepModel, validators
from snappy_pipeline.models.annotation import Vep


class Tool(enum.StrEnum):
    vep = "vep"


class SomaticVariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    path_somatic_variant: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    vep: Vep | None = None
