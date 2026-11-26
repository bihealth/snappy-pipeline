from snappy_pipeline.workflows.adapter_trimming import AdapterTrimmingWorkflow
from snappy_pipeline.workflows.alignment_processing import AlignmentProcessingWorkflow
from snappy_pipeline.workflows.cbioportal_export import cbioportalExportWorkflow
from snappy_pipeline.workflows.combine_variants import CombineVariantsWorkflow
from snappy_pipeline.workflows.concat_variants import ConcatVariantsWorkflow
from snappy_pipeline.workflows.create_proteome import CreateProteomeWorkflow
from snappy_pipeline.workflows.gcnv_model import GcnvModelWorkflow
from snappy_pipeline.workflows.gene_expression_quantification import (
    GeneExpressionQuantificationWorkflow,
)
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from snappy_pipeline.workflows.homologous_recombination_deficiency import (
    HomologousRecombinationDeficiencyWorkflow,
)
from snappy_pipeline.workflows.igv_session_generation import IgvSessionGenerationWorkflow
from snappy_pipeline.workflows.external_data import ExternalDataWorkflow
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
from snappy_pipeline.workflows.somatic_neoepitope_prediction import (
    SomaticNeoepitopePredictionWorkflow,
)
from snappy_pipeline.workflows.somatic_purity_ploidy_estimate import (
    SomaticPurityPloidyEstimateWorkflow,
)
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_signatures import SomaticVariantSignaturesWorkflow
from snappy_pipeline.workflows.somatic_wgs_cnv_calling import SomaticWgsCnvCallingWorkflow
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
from snappy_pipeline.workflows.variant_export_external import VariantExportExternalWorkflow
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow
from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow
from snappy_pipeline.workflows.umi_mapping import UmiMappingWorkflow
from snappy_pipeline.workflows.variant_report import VariantReportWorkflow
from snappy_pipeline.workflows.wgs_cnv_export_external import WgsCnvExportExternalWorkflow
from snappy_pipeline.workflows.wgs_sv_export_external import WgsSvExportExternalWorkflow

#: Step name -> workflow class, from each class's ``name``
WORKFLOW_REGISTRY = {
    cls.name: cls
    for cls in (
        AdapterTrimmingWorkflow,
        AlignmentProcessingWorkflow,
        cbioportalExportWorkflow,
        CombineVariantsWorkflow,
        ConcatVariantsWorkflow,
        CreateProteomeWorkflow,
        GcnvModelWorkflow,
        GeneExpressionQuantificationWorkflow,
        HlaTypingWorkflow,
        HomologousRecombinationDeficiencyWorkflow,
        IgvSessionGenerationWorkflow,
        ExternalDataWorkflow,
        NgsDataQcWorkflow,
        NgsMappingWorkflow,
        PanelOfNormalsWorkflow,
        ReferenceDownloadWorkflow,
        ReferenceIndexWorkflow,
        RepeatExpansionWorkflow,
        SomaticCnvCheckingWorkflow,
        SomaticGeneFusionCallingWorkflow,
        SomaticHlaLohCallingWorkflow,
        SomaticMsiCallingWorkflow,
        SomaticNeoepitopePredictionWorkflow,
        SomaticPurityPloidyEstimateWorkflow,
        SomaticTargetedSeqCnvCallingWorkflow,
        SomaticVariantSignaturesWorkflow,
        SomaticWgsCnvCallingWorkflow,
        SvCallingTargetedWorkflow,
        SvCallingWgsWorkflow,
        MeiWorkflow,
        TumorMutationalBurdenCalculationWorkflow,
        UmiMappingWorkflow,
        VarfishExportWorkflow,
        VariantAnnotationWorkflow,
        VariantCallingWorkflow,
        VariantCheckingWorkflow,
        VariantExportExternalWorkflow,
        VariantFiltrationWorkflow,
        VariantPhasingWorkflow,
        VariantReportWorkflow,
        WgsCnvExportExternalWorkflow,
        WgsSvExportExternalWorkflow,
    )
}
