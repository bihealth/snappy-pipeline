import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel, EnumField


class Tool(enum.StrEnum):
    peddy = "peddy"


class VariantChecking(SnappyStepModel):
    tools_ngs_mapping: list[str] = []
    """copied from ngs mapping config"""

    tools_variant_calling: list[str] = []
    """copied from variant calling config"""

    path_variant_calling: Annotated[str, Field(examples=["../variant_calling"])]
    """Path to variant calling"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.peddy], min_length=1)]
