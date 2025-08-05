import enum

from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel


class FiltrationSchema(enum.StrEnum):
    unfiltered = "unfiltered"
    list = "list"
    sets = "sets"


class FilterSet(enum.StrEnum):
    NO_FILTER = "no_filter"
    DKFZ_ONLY = "dkfz_only"
    DKFZ_AND_EBFILTER = "dkfz_and_ebfilter"
    DKFZ_AND_EBFILTER_AND_OXOG = "dkfz_and_ebfilter_and_oxog"


class TumorMutationalBurden(SnappyStepModel):
    has_annotation: bool = False
    """Needed for building filenames only"""

    path_somatic_variant: Annotated[
        str, Field(examples=["../somatic_variant_annotation", "../somatic_variant_calling"])
    ]
    """Path to variant (directory of vcf files)"""

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    tools_somatic_variant_annotation: list[str] = []
    """default to those configured for somatic_variant_annotation"""

    filtration_schema: FiltrationSchema = FiltrationSchema.list
    """Method of variant filtration (if any)"""

    filter_before_annotation: bool = False
    """Flag if filtration has been done before annotation"""

    # TODO: remove filtration by sets
    # When this is done, the deprecated options can be deleted, and the
    # signature computations will be done just on the output of whatever step
    # is used on input (somatic_variant_calling, somatic_variant_annotation or somatic_variant_filtration)
    filter_sets: Annotated[
        list[FilterSet], Field(None, deprecated="use `filter_list` instead")
    ] = []
    """
    DEPRECATED: use `filter_list` in `somatic_variant_filtration` instead.
    The TMB is computed for all sets decribed here (in conjunction with `exon_lists`).
    When empty, the set of filters is taken from those configured in `somatic_variant_filtration` step.
    """

    exon_lists: Annotated[
        list[str],
        Field(
            "genome_wide",
            deprecated="Works together with filter_set, ignored when `filter_list` is selected",
        ),
    ] = []
    """
    DEPRECATED: use `filter_list` in `somatic_variant_filtration` instead.
    The TMB is computed for all regions decribed here (in conjunction with `filter_sets`).
    When empty, the regions is taken from those configured in `somatic_variant_filtration` step.
    """

    target_regions: str
    """Path to target_regions file (bed format)"""

    missense_regex: str = r".*[\|&]missense_variant[\|&].*"
    """change if the annotation tool doesn't use 'missense_variant' to indicate missense variant"""
