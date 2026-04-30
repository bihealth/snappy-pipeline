from snappy_pipeline.workflows.adapter_trimming import AdapterTrimmingWorkflow
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow

WORKFLOW_REGISTRY = {
    "adapter_trimming": AdapterTrimmingWorkflow,
    "ngs_mapping": NgsMappingWorkflow,
    "somatic_variant_calling": SomaticVariantCallingWorkflow,
}
