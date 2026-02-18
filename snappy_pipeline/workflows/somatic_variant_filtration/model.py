from typing import TypedDict

from pydantic import Field

from snappy_pipeline.models import SnappyModel

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from snappy_pipeline.workflows.any_variant_filtration.model import AnyVariantFiltration
from snappy_pipeline.workflows.any_variant_filtration.model import Filter as ParentFilter


class Ebfilter(SnappyModel):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: int = 20
    min_baseq: int = 15
    path_panel_of_normals_sample_list: str = ""


class EbfilterOnly(TypedDict, total=False):
    ebfilter: Ebfilter


class SomaticVariantFiltration(AnyVariantFiltration):
    filter_list: list[ParentFilter | EbfilterOnly] = []
    variant_origin: VariantOrigin = Field(VariantOrigin.SOMATIC, frozen=True)
