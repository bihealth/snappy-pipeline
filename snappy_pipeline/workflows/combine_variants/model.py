import enum
from typing import Annotated

from pydantic import BaseModel

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema


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
    somatic_variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedSomaticVariant),
    ]
    germline_variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedGermlineVariant),
    ]


class CombineVariants(SnappyStepModel):
    depends_on: CombineVariantsDependsOn
    rename_combined: RenameCombine | None = None
