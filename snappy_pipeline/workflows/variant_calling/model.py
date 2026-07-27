import enum
import logging
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import (
    EnumField,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
)
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.models.parallel import Parallel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

logger = logging.getLogger(__name__)


class ExpectedGermlineVariants(BaseModel):
    """Consumer-driven contract: expected output keys from a variant_calling upstream task."""

    vcf: str
    vcf_tbi: str


class ExpectedSomaticVariants(BaseModel):
    """Consumer-driven contract: expected output keys from a variant_calling upstream task."""

    vcf: str
    vcf_tbi: str


class BafFileGeneration(ToggleModel):
    min_dp: Annotated[int, Field(ge=1)] = 10
    """minimal DP of variant, must be >=1"""


class BcftoolsStats(ToggleModel):
    pass


class JannovarStats(ToggleModel):
    path_ser: str


class BcftoolsRoh(ToggleModel):
    path_targets: str | None = None  # FIXME this says "REQUIRED; optional" in the original code

    path_af_file: str

    ignore_homref: bool = False

    skip_indels: bool = False

    rec_rate: float = 1e-8


class Tool(enum.StrEnum):
    bcftools_call = "bcftools_call"
    gatk3_hc = "gatk3_hc"
    gatk3_ug = "gatk3_ug"
    gatk4_hc_joint = "gatk4_hc_joint"
    gatk4_hc_gvcf = "gatk4_hc_gvcf"
    mutect2 = "mutect2"


class BcftoolsCall(SnappyModel):
    max_depth: int = 250
    max_indel_depth: int = 250
    window_length: int = 10000000
    num_threads: int = 16


class Gatk3Hc(SnappyModel):
    num_threads: int = 16
    window_length: int = 10000000
    allow_seq_dict_incompatibility: bool = False


class Gatk3Ug(SnappyModel):
    num_threads: int = 16
    window_length: int = 10000000
    allow_seq_dict_incompatibility: bool = False
    downsample_to_coverage: int = 250


class Gatk4HcJoint(SnappyModel):
    window_length: int = 10000000
    num_threads: int = 16
    allow_seq_dict_incompatibility: bool = False


class Gatk4HcGvcf(SnappyModel):
    window_length: int = 10000000
    num_threads: int = 16
    allow_seq_dict_incompatibility: bool = False


class TumorNormalMode(enum.StrEnum):
    AUTOMATIC = "automatic"
    PAIRED = "paired"
    TUMOR_ONLY = "tumor_only"
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Contamination(ToggleModel, GATK):
    common_variants: str = ""
    """Common germline variants for contamination estimation"""

    pileup: GATK = GATK()
    """Parameters for GetPileupSummaries used on the tumor bam (& normalk if present)"""

    @model_validator(mode="after")
    def ensure_common_variant_when_enabled(self):
        if self.enabled and not self.common_variants:
            raise ValueError("Common variants must be present when contamination is enabled")
        return self


class Mutect2(Parallel, GATK):
    # Sadly a type of
    # `FilePath | None = None`
    # still applies `FilePath` validation on `None`, which errors
    panel_of_normals: str | None = ""
    """Set path to panel of normals vcf if required"""

    germline_resource: str | None = ""
    """Germline variants resource (same as panel of normals)"""

    contamination: Contamination
    """Estimation of contamination using GetPileupSummaries & CalculateContamination"""

    filtration: GATK = GATK()
    """Additional arguments for filtration"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class VariantCallingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class VariantCalling(SnappyStepModel):
    depends_on: VariantCallingDependsOn = Field(default_factory=VariantCallingDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.gatk4_hc_gvcf)]

    ignore_chroms: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]

    baf_file_generation: BafFileGeneration = BafFileGeneration()

    bcftools_stats: BcftoolsStats | None = None

    jannovar_stats: JannovarStats | None = None

    bcftools_roh: BcftoolsRoh | None = None

    bcftools_call: BcftoolsCall | None = None

    gatk3_hc: Gatk3Hc | None = None

    gatk3_ug: Gatk3Ug | None = None

    gatk4_hc_joint: Gatk4HcJoint | None = None

    gatk4_hc_gvcf: Gatk4HcGvcf | None = None

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""
