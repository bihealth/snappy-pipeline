import enum

from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyStepModel


class SomaticVariantStep(enum.StrEnum):
    CALL = "somatic_variant_calling"
    ANNOTATION = "somatic_variant_annotation"
    FILTER = "somatic_variant_filtration"


class TumorMutationalBurden(SnappyStepModel):
    has_annotation: bool = False
    """Needed for building filenames only"""

    path_somatic_variant: Annotated[
        str, Field(examples=["../somatic_variant_annotation", "../somatic_variant_calling"])
    ]
    """Path to variant (directory of vcf files)"""

    somatic_variant_step: SomaticVariantStep = SomaticVariantStep.FILTER
    """Which pipeline step is used to compute signatures"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    tools_somatic_variant_annotation: list[str] = []
    """default to those configured for somatic_variant_annotation"""

    has_annotation: bool = True
    """Has the inpyut vcf been annotated"""

    is_filtered: bool = True
    """Has the input vcf been filtered"""

    target_regions: str
    """Path to target_regions file (bed format)"""

    missense_regex: str = r".*[\|&]missense_variant[\|&].*"
    """change if the annotation tool doesn't use 'missense_variant' to indicate missense variant"""

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
