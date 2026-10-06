from snappy_pipeline.workflows.any_variant_filtration import AnyVariantFiltrationWorkflow
from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import GermlineVariantFiltration as GermlineVariantFiltrationConfigModel


class GermlineVariantFiltrationWorkflow(AnyVariantFiltrationWorkflow):
    name = "germline_variant_filtration"
    variant_origin = VariantOrigin.GERMLINE
    model_class = GermlineVariantFiltrationConfigModel
