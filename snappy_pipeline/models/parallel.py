import enum

from snappy_pipeline.models import SnappyModel


class Keep(enum.StrEnum):
    ALWAYS = "always"
    NEVER = "never"
    ONERROR = "onerror"


class Parallel(SnappyModel):
    num_cores: int = 2
    """number of cores to use locally"""

    num_jobs: int = 24
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
