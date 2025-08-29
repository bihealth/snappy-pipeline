import enum

from typing import Annotated

from pydantic import Field

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
