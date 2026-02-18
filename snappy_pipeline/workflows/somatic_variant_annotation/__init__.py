from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions

from snappy_pipeline.workflows.any_variant_annotation import AnyVariantAnnotationWorkflow
from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import SomaticVariantAnnotation as SomaticVariantAnnotationConfigModel


class SomaticVariantAnnotationWorkflow(AnyVariantAnnotationWorkflow):
    name = "somatic_variant_annotation"
    variant_origin = VariantOrigin.SOMATIC
    model_class = SomaticVariantAnnotationConfigModel

    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }
