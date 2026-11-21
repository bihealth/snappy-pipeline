import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    RELATED_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedVariants,
)


class Tool(enum.StrEnum):
    peddy = "peddy"


class VariantCheckingDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        RELATED_LEVEL,
        ExpectedPathSchema(ExpectedVariants),
    ]


class VariantChecking(SnappyStepModel, ToolModel):
    depends_on: VariantCheckingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]
