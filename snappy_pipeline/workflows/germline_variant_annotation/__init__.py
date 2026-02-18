from snappy_pipeline.workflows.any_variant_annotation import AnyVariantAnnotationWorkflow
from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import GermlineVariantAnnotation as GermlineVariantAnnotationConfigModel


class GermlineVariantAnnotationWorkflow(AnyVariantAnnotationWorkflow):
    name = "germline_variant_annotation"
    variant_origin = VariantOrigin.GERMLINE
    model_class = GermlineVariantAnnotationConfigModel
