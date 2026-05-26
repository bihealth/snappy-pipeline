import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.somatic_variant_calling.model import ExpectedSomaticVariants


class SomaticVariantStep(enum.StrEnum):
    CALL = "somatic_variant_calling"
    ANNOTATION = "somatic_variant_annotation"
    FILTER = "somatic_variant_filtration"


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

    somatic_variant_step: SomaticVariantStep = SomaticVariantStep.CALL
    """Which pipeline step is used to compute signatures"""

    has_annotation: bool = False
    """Needed for building filenames only"""

    is_filtered: bool = False
    """Has the vcf been already filtered"""

    @model_validator(mode="after")
    def ensure_annotation_and_filtration_are_configured_correctly(self):
        if self.somatic_variant_step == SomaticVariantStep.CALL:
            if self.has_annotation or self.is_filtered:
                raise ValueError(
                    "When the input step is 'somatic_variant_calling', the annotation & filtration status must be set to 'False'"
                )
        elif self.somatic_variant_step == SomaticVariantStep.ANNOTATION:
            if not self.has_annotation:
                raise ValueError(
                    "When the input step is 'somatic_variant_annotation', the annotation status must be set to 'True'"
                )
        elif self.somatic_variant_step == SomaticVariantStep.FILTER:
            if not self.is_filtered:
                raise ValueError(
                    "When the input step is 'somatic_variant_filtration', the filtration status must be set to 'True'"
                )
        return self
