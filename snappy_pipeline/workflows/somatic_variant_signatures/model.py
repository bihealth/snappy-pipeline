import enum

from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel


class SomaticVariantStep(enum.StrEnum):
    CALL = "somatic_variant_calling"
    ANNOTATION = "somatic_variant_annotation"
    FILTER = "somatic_variant_filtration"


class SomaticVariantSignatures(SnappyStepModel):
    path_somatic_variant: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    somatic_variant_step: SomaticVariantStep = SomaticVariantStep.FILTER
    """Which pipeline step is used to compute signatures"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    tools_somatic_variant_annotation: list[str] = []
    """default to those configured for somatic_variant_annotation"""

    has_annotation: bool = False
    """Needed for building filenames only"""

    is_filtered: bool = False
    """Has the vcf been already filtered"""
