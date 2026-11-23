import enum
from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    ANY_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
)


class RenameCombine(enum.StrEnum):
    TUMOR = "tumor"
    GERMLINE = "germline"


class CombineVariantsDependsOn(SnappyModel):
    somatic_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
    ]
    germline_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
        ANY_LEVEL,
    ]
    """The germline calls of the matched normal: its own, or the joint calls of its related group"""

    reference: Reference


class CombineVariants(SnappyStepModel, ToolModel):
    depends_on: CombineVariantsDependsOn

    rename_combined: RenameCombine | None = None
