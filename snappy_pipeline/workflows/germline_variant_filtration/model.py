from pydantic import Field

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from snappy_pipeline.workflows.any_variant_filtration.model import AnyVariantFiltration


class GermlineVariantFiltration(AnyVariantFiltration):
    variant_origin: VariantOrigin = Field(VariantOrigin.GERMLINE, frozen=True)
