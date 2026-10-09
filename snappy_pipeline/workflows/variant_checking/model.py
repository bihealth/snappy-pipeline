import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants


class Tool(enum.StrEnum):
    peddy = "peddy"


class VariantCheckingDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ]


class VariantChecking(SnappyStepModel):
    depends_on: VariantCheckingDependsOn

    """Path to variant calling"""

    tool: Annotated[Tool, EnumField(Tool)]
