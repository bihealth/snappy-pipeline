from snappy_pipeline.workflows.adapter_trimming import AdapterTrimmingWorkflow
from snappy_pipeline.workflows.cbioportal_export import cbioportalExportWorkflow
from snappy_pipeline.workflows.gene_expression_quantification import (
    GeneExpressionQuantificationWorkflow,
)
from snappy_pipeline.workflows.gene_expression_report import GeneExpressionReportWorkflow
from snappy_pipeline.workflows.helper_gcnv_model_targeted import (
    HelperBuildTargetSeqGcnvModelWorkflow,
)
from snappy_pipeline.workflows.helper_gcnv_model_wgs import HelperBuildWgsGcnvModelWorkflow
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from snappy_pipeline.workflows.homologous_recombination_deficiency import (
    HomologousRecombinationDeficiencyWorkflow,
)
from snappy_pipeline.workflows.igv_session_generation import IgvSessionGenerationWorkflow
from snappy_pipeline.workflows.link_in import LinkInWorkflow
from snappy_pipeline.workflows.ngs_data_qc import NgsDataQcWorkflow
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.panel_of_normals import PanelOfNormalsWorkflow
from snappy_pipeline.workflows.reference_download import ReferenceDownloadWorkflow
from snappy_pipeline.workflows.reference_index import ReferenceIndexWorkflow
from snappy_pipeline.workflows.repeat_expansion import RepeatExpansionWorkflow
from snappy_pipeline.workflows.somatic_cnv_checking import SomaticCnvCheckingWorkflow
from snappy_pipeline.workflows.somatic_gene_fusion_calling import SomaticGeneFusionCallingWorkflow
from snappy_pipeline.workflows.somatic_hla_loh_calling import SomaticHlaLohCallingWorkflow
from snappy_pipeline.workflows.somatic_msi_calling import SomaticMsiCallingWorkflow
from snappy_pipeline.workflows.somatic_purity_ploidy_estimate import (
    SomaticPurityPloidyEstimateWorkflow,
)
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow
from snappy_pipeline.workflows.somatic_variant_signatures import SomaticVariantSignaturesWorkflow
from snappy_pipeline.workflows.somatic_wgs_cnv_calling import SomaticWgsCnvCallingWorkflow
from snappy_pipeline.workflows.somatic_wgs_sv_calling import SomaticWgsSvCallingWorkflow
from snappy_pipeline.workflows.sv_calling_targeted import SvCallingTargetedWorkflow
from snappy_pipeline.workflows.sv_calling_wgs import SvCallingWgsWorkflow
from snappy_pipeline.workflows.targeted_seq_mei_calling import MeiWorkflow
from snappy_pipeline.workflows.tumor_mutational_burden import (
    TumorMutationalBurdenCalculationWorkflow,
)
from snappy_pipeline.workflows.varfish_export import VarfishExportWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow
from snappy_pipeline.workflows.variant_checking import VariantCheckingWorkflow
from snappy_pipeline.workflows.variant_denovo_filtration import VariantDeNovoFiltrationWorkflow
from snappy_pipeline.workflows.variant_export_external import VariantExportExternalWorkflow
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow
from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow
from snappy_pipeline.workflows.wgs_cnv_export_external import WgsCnvExportExternalWorkflow
from snappy_pipeline.workflows.wgs_sv_export_external import WgsSvExportExternalWorkflow

WORKFLOW_REGISTRY = {
    "adapter_trimming": AdapterTrimmingWorkflow,
    "cbioportal_export": cbioportalExportWorkflow,
    "gene_expression_quantification": GeneExpressionQuantificationWorkflow,
    "gene_expression_report": GeneExpressionReportWorkflow,
    "helper_gcnv_model_targeted": HelperBuildTargetSeqGcnvModelWorkflow,
    "helper_gcnv_model_wgs": HelperBuildWgsGcnvModelWorkflow,
    "hla_typing": HlaTypingWorkflow,
    "homologous_recombination_deficiency": HomologousRecombinationDeficiencyWorkflow,
    "igv_session_generation": IgvSessionGenerationWorkflow,
    "link_in": LinkInWorkflow,
    "ngs_data_qc": NgsDataQcWorkflow,
    "ngs_mapping": NgsMappingWorkflow,
    "panel_of_normals": PanelOfNormalsWorkflow,
    "reference_download": ReferenceDownloadWorkflow,
    "reference_index": ReferenceIndexWorkflow,
    "repeat_expansion": RepeatExpansionWorkflow,
    "somatic_cnv_checking": SomaticCnvCheckingWorkflow,
    "somatic_gene_fusion_calling": SomaticGeneFusionCallingWorkflow,
    "somatic_hla_loh_calling": SomaticHlaLohCallingWorkflow,
    "somatic_msi_calling": SomaticMsiCallingWorkflow,
    "somatic_purity_ploidy_estimate": SomaticPurityPloidyEstimateWorkflow,
    "somatic_targeted_seq_cnv_calling": SomaticTargetedSeqCnvCallingWorkflow,
    "somatic_variant_calling": SomaticVariantCallingWorkflow,
    "somatic_variant_signatures": SomaticVariantSignaturesWorkflow,
    "somatic_wgs_cnv_calling": SomaticWgsCnvCallingWorkflow,
    "somatic_wgs_sv_calling": SomaticWgsSvCallingWorkflow,
    "sv_calling_targeted": SvCallingTargetedWorkflow,
    "sv_calling_wgs": SvCallingWgsWorkflow,
    "targeted_seq_mei_calling": MeiWorkflow,
    "tumor_mutational_burden": TumorMutationalBurdenCalculationWorkflow,
    "varfish_export": VarfishExportWorkflow,
    "variant_annotation": VariantAnnotationWorkflow,
    "variant_calling": VariantCallingWorkflow,
    "variant_checking": VariantCheckingWorkflow,
    "variant_denovo_filtration": VariantDeNovoFiltrationWorkflow,
    "variant_export_external": VariantExportExternalWorkflow,
    "variant_filtration": VariantFiltrationWorkflow,
    "variant_phasing": VariantPhasingWorkflow,
    "wgs_cnv_export_external": WgsCnvExportExternalWorkflow,
    "wgs_sv_export_external": WgsSvExportExternalWorkflow,
}
