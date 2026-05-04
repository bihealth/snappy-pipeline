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

    """Path to variant calling"""

    tool: Annotated[Tool, EnumField(Tool, default=Tool.peddy)]
