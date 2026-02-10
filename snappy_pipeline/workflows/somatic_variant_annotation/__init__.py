from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions

from snappy_pipeline.workflows.any_variant_annotation import AnyVariantAnnotationWorkflow


class SomaticVariantAnnotationWorkflow(AnyVariantAnnotationWorkflow):
    name = "somatic_variant_annotation"

    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }
