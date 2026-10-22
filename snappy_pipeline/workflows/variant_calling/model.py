import enum
import logging
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import (
    ConfigPath,
    EnumField,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ToolModel,
    validators,
)
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.models.parallel import Parallel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    CommonVariants,
    Dbsnp,
    ExpectedPathSchema,
    GermlineResource,
    Kits,
    Reference,
)
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


class ExpectedPanelOfNormals(BaseModel):
    """Consumer-driven contract: the panel of normals VCF of a mutect2 panel_of_normals task."""

    vcf: str
    vcf_tbi: str


class BafFileGeneration(ToggleModel, ToolModel):
    min_dp: Annotated[int, Field(ge=1)] = 10
    """minimal DP of variant, must be >=1"""


class BcftoolsStats(ToggleModel, ToolModel):
    pass


class BcftoolsRoh(ToggleModel, ToolModel):
    path_af_file: ConfigPath

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


class IgnoreChroms(SnappyModel):
    """The contigs that a caller skips"""

    ignore_chroms: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""


class BcftoolsCall(ToolModel, IgnoreChroms):
    max_depth: int = 250
    max_indel_depth: int = 250
    window_length: int = 10000000
    num_threads: int = 16


class Gatk3Hc(ToolModel, IgnoreChroms):
    num_threads: int = 16
    window_length: int = 10000000
    allow_seq_dict_incompatibility: bool = False


class Gatk3Ug(ToolModel, IgnoreChroms):
    num_threads: int = 16
    window_length: int = 10000000
    allow_seq_dict_incompatibility: bool = False
    downsample_to_coverage: int = 250


class Gatk4HcJoint(ToolModel, IgnoreChroms):
    window_length: int = 10000000
    num_threads: int = 16
    allow_seq_dict_incompatibility: bool = False


class Gatk4HcGvcf(ToolModel, IgnoreChroms):
    window_length: int = 10000000
    num_threads: int = 16
    allow_seq_dict_incompatibility: bool = False


class TumorNormalMode(enum.StrEnum):
    AUTOMATIC = "automatic"
    PAIRED = "paired"
    TUMOR_ONLY = "tumor_only"
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Contamination(ToggleModel, GATK):
    java_options: str = "-Xms4000m -Xmx8000m"
    """Options for the java run-time of ``CalculateContamination``"""

    pileup: GATK = GATK()
    """Parameters for GetPileupSummaries used on the tumor bam (& normalk if present)"""


class Mutect2(Parallel, GATK, IgnoreChroms):
    contamination: Contamination
    """Estimation of contamination using GetPileupSummaries & CalculateContamination"""

    filtration: GATK = GATK()
    """Additional arguments for filtration"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class VariantCallingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    dbsnp: Dbsnp = ""
    """Needed by the GATK callers"""

    panel_of_normals: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"pon", "mutect2"})),
        ExpectedPathSchema(ExpectedPanelOfNormals),
    ] = ""
    """
    Optional for ``tool: mutect2``: a ``panel_of_normals`` task with ``tool: mutect2``,
    whose panel VCF Mutect2 uses to filter recurrent artefacts.
    """

    germline_resource: GermlineResource = ""
    """Optional for ``tool: mutect2``: population allele frequencies (af-only gnomAD VCF)"""

    common_variants: CommonVariants = ""
    """Required for ``tool: mutect2`` with ``contamination.enabled``: common germline variants"""

    kits: Kits = ""
    """Optional: restricts ``bcftools_roh`` to the targets of the (targeted) libraries' kit"""


TOOL_DEPENDENCIES = {
    tool: ("dbsnp",)
    for tool in (Tool.gatk3_hc, Tool.gatk3_ug, Tool.gatk4_hc_joint, Tool.gatk4_hc_gvcf)
}


class VariantCalling(SnappyStepModel):
    depends_on: VariantCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    baf_file_generation: BafFileGeneration = BafFileGeneration()

    bcftools_stats: BcftoolsStats | None = None

    bcftools_roh: BcftoolsRoh | None = None

    bcftools_call: BcftoolsCall | None = None

    gatk3_hc: Gatk3Hc | None = None

    gatk3_ug: Gatk3Ug | None = None

    gatk4_hc_joint: Gatk4HcJoint | None = None

    gatk4_hc_gvcf: Gatk4HcGvcf | None = None

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self

    @model_validator(mode="after")
    def ensure_common_variant_when_enabled(self):
        if self.mutect2 and self.mutect2.contamination.enabled:
            if not self.depends_on.common_variants:
                raise ValueError(
                    "depends_on.common_variants is required when contamination is enabled"
                )
        return self
