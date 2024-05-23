from typing import Annotated

from pydantic import Field

from models import SnappyStepModel


class TumorMutationalBurden(SnappyStepModel):
    has_annotation: bool

    is_filtered: bool

    path_somatic_variant: Annotated[
        str, Field(examples=["../somatic_variant_annotation", "../somatic_variant_calling"])
    ]

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    tools_somatic_variant_annotation: list[str] = []
    """default to those configured for somatic_variant_annotation"""

    filters: list[str] = []
    """
    When using variants after the somatic_variant_filtration step,
    use "no_filter", "dkfz_only", "dkfz_and_ebfilter" or "dkfz_and_ebfilter_and_oxog"
    """

    filtered_regions: list[str] = []
    """When using variants after the somatic_variant_filtration step, use "genome_wide" or "" """

    target_regions: str

    missense_regex: str = r".*[\|&]missense_variant[\|&].*"
    """change if the annotation tool doesn't use 'missense_variant' to indicate missense variant"""
