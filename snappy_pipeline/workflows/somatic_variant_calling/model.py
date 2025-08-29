import enum
from typing import Annotated

from pydantic import AfterValidator, Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    mutect2 = "mutect2"


class Keep(enum.StrEnum):
    ALWAYS = "always"
    NEVER = "never"
    ONERROR = "onerror"


class Parallel(SnappyModel):
    num_cores: int = 2
    """number of cores to use locally"""

    window_length: int = 3500000
    """split input into windows of this size, each triggers a job"""

    padding: int = 5000
    """Padding around scatter-intervals, in bp."""

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


def argument(args: list[str]) -> list[str]:
    def _is_valid_argument(arg: str) -> bool:
        return arg.startswith("--")

    if invalid_args := list(filter(lambda x: not _is_valid_argument(x), args)):
        raise ValueError(f"invalid arguments: {invalid_args}")
    return args


# Adjustment tumor_only mode
class TumorNormalMode(enum.StrEnum):
    AUTOMATIC = "automatic"
    PAIRED = "paired"
    TUMOR_ONLY = "tumor_only"
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Mutect2(Parallel):
    # Sadly a type of
    # `FilePath | None = None`
    # still applies `FilePath` validation on `None`, which errors
    panel_of_normals: str | None = ""
    """Set path to panel of normals vcf if required"""

    germline_resource: str | None = ""
    """Germline variants resource (same as panel of normals)"""

    common_variants: str | None = ""
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

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class SomaticVariantCalling(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [], min_length=1)]
    """List of tools"""

    path_ngs_mapping: str = "../ngs_mapping"
    """Path to ngs_mapping"""

    ignore_chroms: Annotated[
        list[str],
        Field(examples=["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]),
    ] = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]
    """Patterns of contig names to ignore"""

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""
