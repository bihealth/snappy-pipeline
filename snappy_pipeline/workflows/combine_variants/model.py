import enum
from typing import Annotated

from pydantic import BaseModel

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    ANY_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)


class RenameCombine(enum.StrEnum):
    TUMOR = "tumor"
    GERMLINE = "germline"


class ExpectedSomaticVariant(BaseModel):
    vcf: str
    vcf_tbi: str


class ExpectedGermlineVariant(BaseModel):
    vcf: str
    vcf_tbi: str


class CombineVariantsDependsOn(SnappyModel):
    somatic_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedSomaticVariant),
    ]
    germline_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedGermlineVariant),
        ANY_LEVEL,
    ]
    """The germline calls of the matched normal: its own, or the joint calls of its related group"""

    reference: Reference


class CombineVariants(SnappyStepModel, ToolModel):
    depends_on: CombineVariantsDependsOn
    rename_combined: RenameCombine | None = None
