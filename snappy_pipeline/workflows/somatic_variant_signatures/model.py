from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.somatic_variant_calling.model import ExpectedSomaticVariants


class SomaticVariantSignaturesDependsOn(SnappyModel):
    somatic_variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic"})),
        ExpectedPathSchema(ExpectedSomaticVariants),
    ] = "somatic_variant"


class SomaticVariantSignatures(SnappyStepModel):
    depends_on: SomaticVariantSignaturesDependsOn = Field(
        default_factory=SomaticVariantSignaturesDependsOn
    )
