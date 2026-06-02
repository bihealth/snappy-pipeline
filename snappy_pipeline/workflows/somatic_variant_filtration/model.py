from typing import TypedDict

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from snappy_pipeline.workflows.any_variant_filtration.model import AnyVariantFiltration
from snappy_pipeline.workflows.any_variant_filtration.model import Filter as ParentFilter


class Dkfz(SnappyModel):
    pass


class Ebfilter(SnappyModel):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: int = 20
    min_baseq: int = 15
    path_panel_of_normals_sample_list: str = ""


class EbfilterAndDkfz(TypedDict, total=False):
    ebfilter: Ebfilter
    dkfz: Dkfz


class SomaticVariantFiltration(AnyVariantFiltration):
    filter_list: list[ParentFilter | EbfilterAndDkfz] = []
    variant_origin: VariantOrigin = Field(VariantOrigin.SOMATIC, frozen=True)

    @model_validator(mode="after")
    def ensure_filter_list_is_configured_correctly(self):
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
