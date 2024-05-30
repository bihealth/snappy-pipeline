import enum
from typing import Annotated, Any, Self, TypedDict

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class DkfzAndEbfilter(SnappyModel):
    ebfilter_threshold: float = 2.4


class DkfzAndEbfilterAndOxog(SnappyModel):
    vaf_threshold: float = 0.08
    coverage_threshold: float = 5


class DkfzAndOxog(SnappyModel):
    vaf_threshold: float = 0.08
    coverage_threshold: float = 5


class FilterSets(SnappyModel):
    no_filter: str | None = None
    dkfz_only: str | None = None
    dkfz_and_ebfilter: DkfzAndEbfilter | None = None
    dkfz_and_ebfilter_and_oxog: DkfzAndEbfilterAndOxog | None = None
    dkfz_and_oxog: DkfzAndOxog | None = None


class EbfilterSet(SnappyModel):
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: float = 20
    min_baseq: float = 15


class Ebfilter(SnappyModel):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: float = 20
    min_baseq: float = 15


class Dkfz(SnappyModel):
    pass


class Bcftools(SnappyModel):
    include: str = ""
    """Expression to be used in bcftools view --include"""

    exclude: str = ""
    """Expression to be used in bcftools view --exclude"""

    @model_validator(mode="after")
    def ensure_include_or_exclude(self) -> Self:
        if not self.include and not self.exclude:
            raise ValueError("Either include or exclude must be set")
        if self.include and self.exclude:
            raise ValueError("Only one of include or exclude may be set")
        return self

    def keywords(self) -> dict[str, str]:
        if self.include:
            return {"include": self.include}
        elif self.exclude:
            return {"exclude": self.exclude}
        return {}


class Regions(SnappyModel):
    include: str = ""
    """Expression to be used in bcftools view --include"""

    exclude: str = ""
    """Expression to be used in bcftools view --exclude"""

    path_bed: Annotated[str, Field(deprecated="Use `exclude` instead")] = ""
    """Bed file of regions to be considered (variants outside are filtered out)"""

    @model_validator(mode="after")
    def ensure_include_or_exclude(self) -> Self:
        if not any((self.include, self.exclude, self.path_bed)):
            raise ValueError("Either include, exclude or path_bed must be set")
        if sum((bool(self.include), bool(self.exclude), bool(self.path_bed))) > 1:
            raise ValueError("Only one of include, exclude or path_bed may be set")
        return self

    def keywords(self) -> dict[str, str]:
        if self.include:
            return {"include": self.include}
        elif self.exclude:
            return {"exclude": self.exclude}
        elif self.path_bed:
            return {"exclude": self.path_bed}  # path_bed is deprecated and replaced by exclude
        return {}


class Protected(SnappyModel):
    path_bed: str
    """Bed file of regions that should not be filtered out at all."""

    def keywords(self) -> dict[str, str]:
        if self.path_bed:
            return {"path_bed": self.path_bed}
        return {}


class Filter(TypedDict, total=False):
    bcftools: Bcftools
    dkfz: Dkfz
    ebfilter: Ebfilter
    regions: Regions
    protected: Protected


class FiltrationSchema(enum.StrEnum):
    list = "list"
    sets = "sets"


class SomaticVariantFiltration(SnappyStepModel):
    path_somatic_variant: Annotated[
        str, Field(examples=["../somatic_variant_annotation", "../somatic_variant_calling"])
    ] = "../somatic_variant"

    path_ngs_mapping: str = "../ngs_mapping"
    """Needed for dkfz & ebfilter"""

    tools_ngs_mapping: list[str] | None = None
    """Default: use those defined in ngs_mapping step"""

    tools_somatic_variant_calling: list[str] | None = None
    """Default: use those defined in somatic_variant_calling step"""

    tools_somatic_variant_annotation: list[str] | None = None
    """Default: use those defined in somatic_variant_annotation step"""

    has_annotation: bool = True

    filtration_schema: FiltrationSchema = FiltrationSchema.list

    filter_sets: Annotated[FilterSets | None, Field(deprecated="use filter_list instead")] = None

    exon_lists: Annotated[dict[str, Any], Field(deprecated="use filter_list instead")] = {}

    eb_filter: Annotated[
        EbfilterSet | None, Field(deprecated="use filter_list instead")
    ] = EbfilterSet()

    filter_list: list[Filter] = []
    """
    Available filters
    dkfz: {}                                         # Not parametrisable
    ebfilter:
      ebfilter_threshold: 2.4
      shuffle_seed: 1
      panel_of_normals_size: 25
      min_mapq: 20
      min_baseq: 15
    bcftools:
      include: ""                                   # Expression to be used in bcftools view --include
      exclude: ""                                   # Expression to be used in bcftools view --exclude
    regions:
      path_bed: REQUIRED                            # Bed file of regions to be considered (variants outside are filtered out)
    protected:
      path_bed: REQUIRED                            # Bed file of regions that should not be filtered out at all.
    """

    @model_validator(mode="after")
    def ensure_filter_list_is_configured_correctly(self):
        if self.filter_list:
            # check ebfilter and dkfz are only used at most once
            num_ebfilter = num_dkfz = 0
            for f in self.filter_list:
                if "ebfilter" in f:
                    num_ebfilter += 1
                if "dkfz" in f:
                    num_dkfz += 1
            if num_ebfilter > 1:
                raise ValueError("Only one ebfilter is allowed")
            if num_dkfz > 1:
                raise ValueError("Only one dkfz is allowed")
        return self

    @model_validator(mode="after")
    def ensure_either_filter_sets_or_filter_list_is_configured(self):
        if self.filtration_schema == FiltrationSchema.sets:
            if not self.filter_sets:
                raise ValueError("filter_sets must be set")
        if self.filtration_schema == FiltrationSchema.list:
            if not self.filter_list:
                raise ValueError("filter_list must be set")
        if self.filter_sets and self.filter_list:
            raise ValueError("Either filter_sets or filter_list must be set")
        return self
