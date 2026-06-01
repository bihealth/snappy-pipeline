import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants


class Tool(enum.StrEnum):
    peddy = "peddy"


class VariantCheckingDependsOn(SnappyModel):
    variant_calling: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ] = "variant_calling"


class VariantChecking(SnappyStepModel):
    depends_on: VariantCheckingDependsOn = Field(default_factory=VariantCheckingDependsOn)

    """Path to variant calling"""

    tool: Annotated[Tool, EnumField(Tool, default=Tool.peddy)]
