import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class ExpectedVariantVcf(SnappyModel):
    """Generic VCF contract consumed by the unified variant_annotation step."""

    vcf: str
    vcf_tbi: str


class ExpectedAnnotatedVariants(SnappyModel):
    """Consumer-driven contract: expected annotated-variant output keys."""

    vcf: str
    vcf_tbi: str


class VariantAnnotationDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ]

    reference: Reference


class VariantAnnotation(SnappyStepModel):
    depends_on: VariantAnnotationDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    vep: Vep = Field(default_factory=Vep)

    mehari: Mehari | None = None
