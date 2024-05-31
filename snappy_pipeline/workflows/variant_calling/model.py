import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToggleModel, validators


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


class VariantCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.gatk4_hc_gvcf], min_length=1)]

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
