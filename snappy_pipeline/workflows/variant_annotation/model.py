import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, LibrarySelectionMixin, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema


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
    variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ]


class VariantAnnotation(LibrarySelectionMixin, SnappyStepModel):
    depends_on: VariantAnnotationDependsOn

    tool: Annotated[Tool, EnumField(Tool, default=Tool.vep)]

    vep: Vep = Field(default_factory=Vep)

    mehari: Mehari | None = None
