import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyStepModel, validators
from snappy_pipeline.models.annotation import Vep


class Tool(enum.StrEnum):
    vep = "vep"


class FiltrationSchema(enum.StrEnum):
    unfiltered = "unfiltered"
    list = "list"
    sets = "sets"


class FilterSets(enum.StrEnum):
    NO_FILTER = "no_filter"
    DKFZ_ONLY = "dkfz_only"
    DKFZ_AND_EBFILTER = "dkfz_and_ebfilter"
    DKFZ_AND_EBFILTER_AND_OXOG = "dkfz_and_ebfilter_and_oxog"


class SomaticVariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    path_somatic_variant: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    filtration_schema: FiltrationSchema = FiltrationSchema.unfiltered
    """Method of variant filtration (if any)"""

    filter_sets: Annotated[FilterSets | None, Field(deprecated="use filter_list instead")] = None

    exon_lists: Annotated[str, Field(deprecated="use filter_list instead")] | None = None

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    vep: Vep | None = None

    @model_validator(mode="after")
    def ensure_filter_sets_properly_configured(self):
        if self.filtration_schema == FiltrationSchema.sets:
            if not self.filter_sets:
                raise ValueError(
                    "filter_sets and exon_lists must be set (set exon_lists to 'genome_wide' if no exon list is provided)"
                )
        return self
