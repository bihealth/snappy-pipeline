import enum
import typing
from typing import Annotated
from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.cnvkit import Cnvkit as CnvkitGeneric
from snappy_pipeline.models.purecn import IntervalFilter
from snappy_pipeline.models.purecn import Segmentation as PureCNSegmentation
from snappy_pipeline.models.purecn import PureCN as PureCNBase
from snappy_pipeline.models.purecn import Variant as PureCNVariantParams
from snappy_pipeline.models.common import LibraryKitEntry, Sex


class WgsCaller(enum.StrEnum):
    CNVKIT = "cnvkit"
    CONTROL_FREEC = "control_freec"


class WesCaller(enum.StrEnum):
    CNVKIT = "cnvkit"
    PURECN = "purecn"
    SEQUENZA = "sequenza"


class Tools(SnappyModel):
    wgs: Annotated[list[WgsCaller], EnumField(WgsCaller, [])]
    """WGS calling tools"""

    wes: Annotated[list[WesCaller], EnumField(WesCaller, [])]
    """WES calling tools"""


class SequencingMethod(enum.StrEnum):
    WES = "hybrid"
    PANEL = "amplicon"
    WGS = "wgs"


class PanelOfNormalsOrigin(enum.StrEnum):
    COHORT = "cohort"
    """Use (& possibly create) a panel of normals from the current cohort of normals in the panel_of_normals step"""

    FILE = "file"
    """Use an panel of normals from another cohort or from public data"""

    FLAT = "flat"
    """Use a flat panel of normal (no panel of normals, actually)"""

    PAIRED = "paired"
    """Use the paired normal as reference (no panel of normal, actually)"""


class PanelOfNormals(SnappyModel):
    enabled: bool = False
    """Use panel of normals during CNV calling"""

    source: PanelOfNormalsOrigin = PanelOfNormalsOrigin.COHORT
    """Which type of panel of normals should be used, cohort is generally recommended"""

    path_panel_of_normals: str = "../panel_of_normals"
    """Path to panel of normals (used for cohort & file sources)"""

    @model_validator(mode="after")
    def ensure_panel_of_normals_path(self):
        if (
            self.enabled
            and self.source != PanelOfNormalsOrigin.FLAT
            and not self.path_panel_of_normals
        ):
            raise ValueError("Undefined panel of normal path")
        return self


class VariantOrigin(enum.StrEnum):
    COHORT = "cohort"
    """Call somatic variants from the current cohort of normals in the somatic_variant_calling step"""

    FILE = "file"
    """Use an panel of normals from another cohort or from public data"""


class VariantTool(enum.StrEnum):
    MUTECT2 = "mutect2"


class Variant(SnappyModel):
    enabled: bool = False
    """Use variants (somatic &/or germline) to improve CNV calling"""

    source: VariantOrigin = VariantOrigin.FILE
    """Where are the variants obtained from"""

    path_somatic_variant_calling: Annotated[
        str,
        Field(
            examples=[
                "../somatic_variant_calling",
                "../somatic_variant_calling_for_CNV",
                "/public_data/common_variants.vcf.gz",
            ]
        ),
    ] = ""
    """
    Path to the variants to use for CNV calling.

    The path can be either to the ``somatic_variant_calling`` step in the pipeline, if "cohort" is selected,
    or to the vcf file with the variants when "file" is selected as source.
    """

    tool: VariantTool = VariantTool.MUTECT2
    """Tool used to call somatic variants in the pipeline"""

    @model_validator(mode="after")
    def ensure_path_to_variants(self):
        if (
            self.enabled
            and self.source == VariantOrigin.FILE
            and not self.path_somatic_variant_calling
        ):
            raise ValueError(
                "A path to the variant vcf file must be provided when selecting 'file' as source"
            )
        return self


class Sequenza(SnappyModel):
    pass


class ControlFreec(SnappyModel):
    pass


class VariantPureCN(Variant, PureCNVariantParams):
    pass


class PureCN(PureCNBase):
    panel_of_normals: PanelOfNormals = PanelOfNormals()
    """
    Panel of normals created by the NormalDB.R script.
    This is required even if the normal/tumor paired mode won't use it.
    """

    @model_validator(mode="after")
    def restrict_pon_mode(self) -> typing.Self:
        if not self.panel_of_normals.enabled:
            raise ValueError("PureCN requires a panel of normals")
        return self

    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect)"""

    sample_sex: Sex = Sex()

    somatic_variant_calling: VariantPureCN = VariantPureCN()

    path_container: Annotated[
        str, Field(examples=["../panel_of_normals/work/containers/out/purecn.simg"])
    ] = ""
    """Conda installation not working well, container is required. When missing the container is downloaded"""


class PurityTool(enum.StrEnum):
    ASCAT = "ascat"
    PURECN = "purecn"


class PurityOrigin(enum.StrEnum):
    AUTOMATIC = "auto"
    """Use current tool to compute purity & ploidy (PureCn & squenza estimate purity & ploidy)"""

    COHORT = "cohort"
    """Use external tool from the pipleine to compute purity & ploidy"""

    SAMPLESHEET = "samplesheet"
    """Extract purity/ploidy from sample sheet"""

    CONFIG = "config"
    """Extract purity/ploidy from configuration file (all samples have the same value)"""


class Purity(SnappyModel):
    enabled: bool = False
    """Use sample purity during CNV calling"""

    source: PurityOrigin = PurityOrigin.SAMPLESHEET

    path_somatic_purity_ploidy_estimate: str = "../somatic_purity_ploidy_estimate"

    tool: PurityTool = PurityTool.PURECN
    """Tool used for purity estimation, if not set, try samplesheet, otherwise default_value"""

    purity: float | None = None
    """Default purity estimate"""
    ploidy: float = 2.0
    """Default ploidy value"""

    @model_validator(mode="after")
    def ensure_valid_params_for_source(self):
        if self.enabled and self.source == PurityOrigin.CONFIG and self.purity is None:
            raise ValueError("Missing default purity value")
        return self


class PanelOfNormalsCnvkit(PanelOfNormals):
    enabled: bool = True
    """Reset enabled value, cnvkit always needs a panel of normal (even flat or paired)"""
    path_targets: str | None = None
    """Path to target file (used only when pon is obtained from file, taken from pipeline step otherwise)"""
    path_antitargets: str | None = None
    """Path to antitarget file (used only when pon is obtained from file, taken from pipeline step otherwise)"""

    @model_validator(mode="after")
    def ensure_paths_target_antitarget(self):
        if self.source == PanelOfNormalsOrigin.FILE:
            if self.path_targets is None or self.path_antitargets is None:
                raise ValueError(
                    "When using a previous pon, target & antitarget files must be defined"
                )
        return self


class VariantCnvkit(Variant):
    min_variant_depth: int = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency plot."""
    zygocity_freq: float | None = None
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""


class Cnvkit(CnvkitGeneric):
    panel_of_normals: PanelOfNormalsCnvkit = PanelOfNormalsCnvkit()

    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """
    Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect)
    """

    somatic_variant_calling: VariantCnvkit = VariantCnvkit()

    somatic_purity_ploidy_estimate: Purity = Purity()

    @model_validator(mode="after")
    def ensure_purity_not_auto(self):
        if self.somatic_purity_ploidy_estimate.source == PurityOrigin.AUTOMATIC:
            raise ValueError("Cnvkit cannot compute purity/ploidy by itself")
        return self

    sample_sex: Sex = Sex()

    path_access: str | None = None
    """Overrides access when not None"""


class SomaticCnvCalling(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"
    """Path to bam files"""

    tools: Tools
    """Tools for WGS & WES data"""

    cnvkit: Cnvkit | None = None
    purecn: PureCN | None = None
    sequenza: Sequenza | None = None
    control_freec: ControlFreec | None = None
