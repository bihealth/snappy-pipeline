import enum
from os import PathLike
from typing import Annotated, Self

from pydantic import AfterValidator, DirectoryPath, Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.Enum):
    MUTECT = "mutect"
    MUTECT2 = "mutect2"
    SCALPEL = "scalpel"


class BcfToolsJoint(SnappyModel):
    max_depth: int = 4000
    max_indel_depth: int = 4000
    window_length: int = 10000000
    num_threads: int = 16


class PlatypusJoint(SnappyModel):
    split_complex_mnvs: bool = True
    """whether or not to split complex and MNV variants"""

    num_threads: int = 16


class Keep(enum.Enum):
    ALWAYS = "always"
    NEVER = "never"
    ONERROR = "onerror"


class Parallel(SnappyModel):
    num_cores: int = 2
    """number of cores to use locally"""

    window_length: int = 3500000
    """split input into windows of this size, each triggers a job"""

    num_jobs: int = 500
    """number of windows to process in parallel"""

    use_profile: bool = True
    """use Snakemake profile for parallel processing"""

    restart_times: int = 5
    """number of times to re-launch jobs in case of failure"""

    max_jobs_per_second: int = 2
    """throttling of job creation"""

    max_status_checks_per_second: int = 10
    """throttling of status checks"""

    debug_trunc_tokens: int = 0
    """truncation to first N tokens (0 for none)"""

    keep_tmpdir: Keep = Keep.NEVER
    """keep temporary directory, {always, never, onerror}"""

    job_mult_memory: float = 1
    """memory multiplier"""

    job_mult_time: float = 1
    """running time multiplier"""

    merge_mult_memory: float = 1
    """memory multiplier for merging"""

    merge_mult_time: float = 1
    """running time multiplier for merging"""


class Mutect(Parallel):
    pass


def argument(args: list[str]) -> list[str]:
    def _is_valid_argument(arg: str) -> bool:
        return arg.startswith("--")

    if any(_is_valid_argument(arg) for arg in args):
        raise ValueError(
            f"invalid arguments: {list(filter(lambda x: not _is_valid_argument(x), args))}"
        )
    return args


class Mutect2(Parallel):

    # Sadly a type of
    # `FilePath | None = None`
    # still applies `FilePath` validation on `None`, which errors
    panel_of_normals: PathLike | None = None
    """Set path to panel of normals vcf if required"""

    germline_resource: PathLike | None = None
    """Germline variants resource (same as panel of normals)"""

    common_variants: PathLike | None = None
    """Common germline variants for contamination estimation"""

    extra_arguments: Annotated[
        list[str],
        AfterValidator(argument),
        Field(
            examples=[
                "--read-filter CigarContainsNoNOperator",
                "--annotation AssemblyComplexity BaseQuality",
            ]
        ),
    ] = []
    """
    List additional Mutect2 arguments.
    Each additional argument must be of the form:
    "--<argument name> <argument value>"
    For example, to filter reads prior to calling & to add annotations to the output vcf:
      - "--read-filter CigarContainsNoNOperator"
      - "--annotation AssemblyComplexity BaseQuality"
    """

    window_length: int = 50000000


class Scalpel(SnappyModel):
    path_target_regions: PathLike


class Strelka2(SnappyModel):
    path_target_regions: PathLike | None = None
    """For exomes: include a bgzipped bed file with tabix index. That also triggers the --exome flag"""


class GatkHcJoint(Parallel):
    # GATK HC--specific configuration
    allow_seq_dict_incompatibility: bool = False
    annotations: list[str] = [
        "BaseQualityRankSumTest",
        "FisherStrand",
        "GCContent",
        "HaplotypeScore",
        "HomopolymerRun",
        "MappingQualityRankSumTest",
        "MappingQualityZero",
        "QualByDepth",
        "ReadPosRankSumTest",
        "RMSMappingQuality",
        "DepthPerAlleleBySample",
        "Coverage",
        "ClippingRankSumTest",
        "DepthPerSampleHC",
    ]

    window_length: int = 50000000


class GatkUgJoint(Parallel):
    # GATK UG--specific configuration
    downsample_to_coverage: float = 250
    allow_seq_dict_incompatibility: bool = False
    annotations: list[str] = [
        "BaseQualityRankSumTest",
        "FisherStrand",
        "GCContent",
        "HaplotypeScore",
        "HomopolymerRun",
        "MappingQualityRankSumTest",
        "MappingQualityZero",
        "QualByDepth",
        "ReadPosRankSumTest",
        "RMSMappingQuality",
        "DepthPerAlleleBySample",
        "Coverage",
        "ClippingRankSumTest",
        "DepthPerSampleHC",
    ]

    window_length: int = 50000000


class SamtoolsMpileup(SnappyModel):
    max_depth: int = 4000

    max_indel_depth: int = 4000

    min_bq: int = 13

    no_baq: bool = True


class VarscanJoint(Parallel, SamtoolsMpileup):
    min_coverage: float = 8

    min_reads2: int = 2

    min_avg_qual: float = 15

    min_var_freq: Annotated[float, Field(ge=0, le=1)] = 0.01

    min_freq_for_hom: Annotated[float, Field(ge=0, le=1)] = 0.75

    p_value: Annotated[float, Field(ge=0, le=1)] = 99e-02

    window_length: int = 5000000


class SomaticVariantCalling(SnappyStepModel):
    tools: Annotated[list[Tool], EnumField(Tool, [])]
    """List of tools"""

    path_ngs_mapping: DirectoryPath
    """Path to ngs_mapping"""

    ignore_chroms: Annotated[
        list[str],
        Field(examples=["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]),
    ] = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]
    """Patterns of contig names to ignore"""

    bcftools_joint: BcfToolsJoint | None = None
    """Configuration for joint calling with samtools+bcftools."""

    platypus_joint: PlatypusJoint | None = None
    """Configuration for joint calling with Platypus."""

    mutect: Mutect | None = None
    """Configuration for MuTect"""

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""

    scalpel: Scalpel | None = None
    """Configuration for Scalpel"""

    strelka2: Strelka2 | None = None
    """Configuration for Strelka2"""

    gatk_hc_joint: GatkHcJoint | None = None
    """Configuration for GatkHcJoint"""

    gatk_ug_joint: GatkUgJoint | None = None
    """Configuration for GatkUgJoint"""

    varscan_joint: VarscanJoint | None = None
    """Configuration for VarscanJoint"""

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
