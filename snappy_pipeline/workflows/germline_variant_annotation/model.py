from pydantic import Field

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from snappy_pipeline.workflows.any_variant_annotation.model import AnyVariantAnnotation


class GermlineVariantAnnotation(AnyVariantAnnotation):
    variant_origin: VariantOrigin = Field(VariantOrigin.GERMLINE, frozen=True)
