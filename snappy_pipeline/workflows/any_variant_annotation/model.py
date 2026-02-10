import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyStepModel, validators
from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from snappy_pipeline.models.annotation import Vep


class Tool(enum.StrEnum):
    vep = "vep"


class AnyVariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    variant_origin: VariantOrigin = VariantOrigin.SOMATIC

    path_variant: Annotated[
        str, Field(examples=["../somatic_variant_calling"], alias="path_somatic_variant")
    ]

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_variant_calling: Annotated[list[str], Field(alias="tools_somatic_variant_calling")] = []
    """default to those configured for somatic_variant_calling"""

    vep: Vep | None = None

    @model_validator(mode="after")
    def ensure_tools_defined_for_any_origin(self):
        if self.variant_origin == VariantOrigin.ANY and not self.tools_variant_calling:
            raise ValueError("Variant calling tools must be defined for 'any' variant origin")
        return self
