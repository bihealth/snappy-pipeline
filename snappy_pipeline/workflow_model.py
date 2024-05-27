import enum
from os import PathLike
from typing import TypedDict

from pydantic import ConfigDict

from snappy_pipeline.models import SnappyModel
from snappy_pipeline.workflows.adapter_trimming.model import AdapterTrimming
from snappy_pipeline.workflows.cbioportal_export.model import CbioportalExport
from snappy_pipeline.workflows.gene_expression_quantification.model import (
    GeneExpressionQuantification,
)
from snappy_pipeline.workflows.gene_expression_report.model import GeneExpressionReport
from snappy_pipeline.workflows.helper_gcnv_model_targeted.model import HelperGcnvModelTargeted
from snappy_pipeline.workflows.helper_gcnv_model_wgs.model import HelperGcnvModelWgs
from snappy_pipeline.workflows.hla_typing.model import HlaTyping
from snappy_pipeline.workflows.homologous_recombination_deficiency.model import (
    HomologousRecombinationDeficiency,
)
from snappy_pipeline.workflows.igv_session_generation.model import IgvSessionGeneration
from snappy_pipeline.workflows.ngs_data_qc.model import NgsDataQc
from snappy_pipeline.workflows.ngs_mapping.model import NgsMapping
from snappy_pipeline.workflows.panel_of_normals.model import PanelOfNormals
from snappy_pipeline.workflows.repeat_expansion.model import RepeatExpansion
from snappy_pipeline.workflows.somatic_cnv_checking.model import SomaticCnvChecking
from snappy_pipeline.workflows.somatic_gene_fusion_calling.model import SomaticGeneFusionCalling
from snappy_pipeline.workflows.somatic_hla_loh_calling.model import SomaticHlaLohCalling
from snappy_pipeline.workflows.somatic_msi_calling.model import SomaticMsiCalling
from snappy_pipeline.workflows.somatic_purity_ploidy_estimate.model import (
    SomaticPurityPloidyEstimate,
)
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling.model import (
    SomaticTargetedSeqCnvCalling,
)
from snappy_pipeline.workflows.somatic_variant_annotation.model import SomaticVariantAnnotation
from snappy_pipeline.workflows.somatic_variant_calling.model import SomaticVariantCalling
from snappy_pipeline.workflows.somatic_variant_filtration.model import SomaticVariantFiltration
from snappy_pipeline.workflows.somatic_variant_signatures.model import SomaticVariantSignatures
from snappy_pipeline.workflows.somatic_wgs_cnv_calling.model import SomaticWgsCnvCalling
from snappy_pipeline.workflows.somatic_wgs_sv_calling.model import SomaticWgsSvCalling
from snappy_pipeline.workflows.sv_calling_targeted.model import SvCallingTargeted
from snappy_pipeline.workflows.sv_calling_wgs.model import SvCallingWgs
from snappy_pipeline.workflows.targeted_seq_mei_calling.model import TargetedSeqMeiCalling
from snappy_pipeline.workflows.tumor_mutational_burden.model import TumorMutationalBurden
from snappy_pipeline.workflows.varfish_export.model import VarfishExport
from snappy_pipeline.workflows.variant_annotation.model import VariantAnnotation
from snappy_pipeline.workflows.variant_calling.model import VariantCalling
from snappy_pipeline.workflows.variant_checking.model import VariantChecking
from snappy_pipeline.workflows.variant_denovo_filtration.model import VariantDenovoFiltration
from snappy_pipeline.workflows.variant_export_external.model import VariantExportExternal
from snappy_pipeline.workflows.variant_filtration.model import VariantFiltration
from snappy_pipeline.workflows.variant_phasing.model import VariantPhasing
from snappy_pipeline.workflows.wgs_cnv_export_external.model import WgsCnvExportExternal
from snappy_pipeline.workflows.wgs_sv_export_external.model import WgsSvExportExternal


class PathModel(SnappyModel):
    path: PathLike = ""


class StaticDataConfig(SnappyModel):
    reference: PathModel
    cosmic: PathModel | None
    dbsnp: PathModel | None


class SearchPattern(SnappyModel):
    left: str = "*.R1.fastq.gz"
    right: str | None = "*.R2.fastq.gz"


class DataSetType(enum.StrEnum):
    MATCHED_CANCER = "matched_cancer"
    GERMLINE_VARIANTS = "germline_variants"


class NamingScheme(enum.StrEnum):
    ONLY_SECONDARY_ID = "only_secondary_id"


class DataSet(SnappyModel):
    file: PathLike = ""
    search_patterns: list[SearchPattern] = [SearchPattern()]
    search_paths: list[PathLike] = ["../raw"]
    type: DataSetType = DataSetType.MATCHED_CANCER
    naming_scheme: NamingScheme = NamingScheme.ONLY_SECONDARY_ID


class StepConfig(TypedDict, total=False):
    adapter_trimming: AdapterTrimming
    cbioportal_export: CbioportalExport
    gene_expression_quantification: GeneExpressionQuantification
    gene_expression_report: GeneExpressionReport
    helper_gcnv_model_targeted: HelperGcnvModelTargeted
    helper_gcnv_model_wgs: HelperGcnvModelWgs
    hla_typing: HlaTyping
    homologous_recombination_deficiency: HomologousRecombinationDeficiency
    igv_session_generation: IgvSessionGeneration
    ngs_data_qc: NgsDataQc
    ngs_mapping: NgsMapping
    panel_of_normals: PanelOfNormals
    repeat_expansion: RepeatExpansion
    somatic_cnv_checking: SomaticCnvChecking
    somatic_gene_fusion_calling: SomaticGeneFusionCalling
    somatic_hla_loh_calling: SomaticHlaLohCalling
    somatic_msi_calling: SomaticMsiCalling
    somatic_purity_ploidy_estimate: SomaticPurityPloidyEstimate
    somatic_targeted_seq_cnv: SomaticTargetedSeqCnvCalling
    somatic_variant_annotation: SomaticVariantAnnotation
    somatic_variant_calling: SomaticVariantCalling
    somatic_variant_filtration: SomaticVariantFiltration
    somatic_variant_signatures: SomaticVariantSignatures
    somatic_wgs_cnv_calling: SomaticWgsCnvCalling
    somatic_wgs_sv_calling: SomaticWgsSvCalling
    sv_calling_targeted: SvCallingTargeted
    sv_calling_wgs: SvCallingWgs
    targeted_seq_mei_calling: TargetedSeqMeiCalling
    tumor_mutational_burden: TumorMutationalBurden
    varfish_export: VarfishExport
    variant_annotation: VariantAnnotation
    variant_calling: VariantCalling
    variant_checking: VariantChecking
    variant_denovo_filtration: VariantDenovoFiltration
    variant_export_external: VariantExportExternal
    variant_filtration: VariantFiltration
    variant_phasing: VariantPhasing
    wgs_cnv_export: WgsCnvExportExternal
    wgs_sv_export: WgsSvExportExternal


class ConfigModel(SnappyModel):
    model_config = ConfigDict(
        extra="allow",
        use_attribute_docstrings=True,
        use_enum_values=True,
    )

    static_data_config: StaticDataConfig
    step_config: StepConfig
    data_sets: dict[str, DataSet]
