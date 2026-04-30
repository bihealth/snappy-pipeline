import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    peddy = "peddy"


class VariantCheckingDependsOn(SnappyModel):
    variant_calling: str = "variant_calling"


class VariantChecking(SnappyStepModel):
    depends_on: VariantCheckingDependsOn = Field(default_factory=VariantCheckingDependsOn)

    tools_ngs_mapping: list[str] = []
    """copied from ngs mapping config"""

    tools_variant_calling: list[str] = []
    """copied from variant calling config"""

    """Path to variant calling"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.peddy], min_length=1)]
