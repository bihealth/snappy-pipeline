"""The workflow steps: step name -> workflow class.

A step's module is ``snappy_pipeline.workflows.<step name>``.
The class is imported on first access, so listing the step names (the CLI and its shell completion)
imports no step module.
"""

import importlib
from collections.abc import Iterator, Mapping, MutableMapping
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snappy_pipeline.workflows.abstract import BaseStep

#: Step name -> the name of its workflow class in ``snappy_pipeline.workflows.<step name>``
WORKFLOW_CLASSES = {
    "adapter_trimming": "AdapterTrimmingWorkflow",
    "alignment_processing": "AlignmentProcessingWorkflow",
    "cbioportal_export": "cbioportalExportWorkflow",
    "combine_variants": "CombineVariantsWorkflow",
    "concat_variants": "ConcatVariantsWorkflow",
    "create_proteome": "CreateProteomeWorkflow",
    "gcnv_model": "GcnvModelWorkflow",
    "gene_expression_quantification": "GeneExpressionQuantificationWorkflow",
    "hla_typing": "HlaTypingWorkflow",
    "homologous_recombination_deficiency": "HomologousRecombinationDeficiencyWorkflow",
    "igv_session_generation": "IgvSessionGenerationWorkflow",
    "external_data": "ExternalDataWorkflow",
    "ngs_data_qc": "NgsDataQcWorkflow",
    "ngs_mapping": "NgsMappingWorkflow",
    "panel_of_normals": "PanelOfNormalsWorkflow",
    "reference_download": "ReferenceDownloadWorkflow",
    "reference_index": "ReferenceIndexWorkflow",
    "repeat_expansion": "RepeatExpansionWorkflow",
    "somatic_cnv_checking": "SomaticCnvCheckingWorkflow",
    "somatic_gene_fusion_calling": "SomaticGeneFusionCallingWorkflow",
    "somatic_hla_loh_calling": "SomaticHlaLohCallingWorkflow",
    "somatic_msi_calling": "SomaticMsiCallingWorkflow",
    "somatic_neoepitope_prediction": "SomaticNeoepitopePredictionWorkflow",
    "somatic_purity_ploidy_estimate": "SomaticPurityPloidyEstimateWorkflow",
    "somatic_targeted_seq_cnv_calling": "SomaticTargetedSeqCnvCallingWorkflow",
    "somatic_variant_signatures": "SomaticVariantSignaturesWorkflow",
    "somatic_wgs_cnv_calling": "SomaticWgsCnvCallingWorkflow",
    "sv_calling_targeted": "SvCallingTargetedWorkflow",
    "sv_calling_wgs": "SvCallingWgsWorkflow",
    "targeted_seq_mei_calling": "MeiWorkflow",
    "tumor_mutational_burden": "TumorMutationalBurdenCalculationWorkflow",
    "umi_mapping": "UmiMappingWorkflow",
    "varfish_export": "VarfishExportWorkflow",
    "variant_annotation": "VariantAnnotationWorkflow",
    "variant_calling": "VariantCallingWorkflow",
    "variant_checking": "VariantCheckingWorkflow",
    "variant_export_external": "VariantExportExternalWorkflow",
    "variant_filtration": "VariantFiltrationWorkflow",
    "variant_phasing": "VariantPhasingWorkflow",
    "variant_report": "VariantReportWorkflow",
    "wgs_cnv_export_external": "WgsCnvExportExternalWorkflow",
    "wgs_sv_export_external": "WgsSvExportExternalWorkflow",
}


class WorkflowRegistry(MutableMapping[str, type["BaseStep"]]):
    """Step name -> workflow class, imported from the step's module on first access.

    Tests put stub classes in with ``monkeypatch.setitem``.
    """

    def __init__(self, class_names: Mapping[str, str]):
        self._entries: dict[str, str | type[BaseStep]] = dict(class_names)

    def __getitem__(self, step: str) -> type[BaseStep]:
        entry = self._entries[step]
        if isinstance(entry, str):
            module = importlib.import_module(f"snappy_pipeline.workflows.{step}")
            entry = self._entries[step] = getattr(module, entry)
        return entry

    def __setitem__(self, step: str, workflow_cls: type[BaseStep]) -> None:
        self._entries[step] = workflow_cls

    def __delitem__(self, step: str) -> None:
        del self._entries[step]

    def __iter__(self) -> Iterator[str]:
        return iter(self._entries)

    def __len__(self) -> int:
        return len(self._entries)

    def __contains__(self, step: object) -> bool:
        return step in self._entries


#: Step name -> workflow class
WORKFLOW_REGISTRY = WorkflowRegistry(WORKFLOW_CLASSES)
