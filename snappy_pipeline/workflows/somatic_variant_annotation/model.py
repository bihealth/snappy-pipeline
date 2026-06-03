import enum
from typing import Annotated

from pydantic import AliasChoices, BaseModel, Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class ExpectedVariantVcf(BaseModel):
    vcf: str
    vcf_tbi: str


class SomaticVariantAnnotationDependsOn(SnappyModel):
    variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ] = Field(default="", validation_alias=AliasChoices("variant", "somatic_variant"))


class SomaticVariantAnnotation(SnappyStepModel):
    depends_on: SomaticVariantAnnotationDependsOn = Field(
        default_factory=SomaticVariantAnnotationDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.vep)]

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    vep: Vep | None = None

    mehari: Mehari | None = None
