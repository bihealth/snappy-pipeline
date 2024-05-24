from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel, SnappyModel, KeepTmpdir


class GatkReadBackedPhasing(SnappyModel):
    phase_quality_threshold: float = 20.0
    """quality threshold for phasing"""

    window_length: int = 5000000
    """split input into windows of this size, each triggers a job"""

    num_jobs: int = 1000
    """number of windows to process in parallel"""

    use_profile: bool = True
    """use Snakemake profile for parallel processing"""

    restart_times: int = 0
    """number of times to re-launch jobs in case of failure"""

    max_jobs_per_second: int = 10
    """throttling of job creation"""

    max_status_checks_per_second: int = 10
    """throttling of status checks"""

    debug_trunc_tokens: int = 0
    """truncation to first N tokens (0 for none)"""

    keep_tmpdir: KeepTmpdir = KeepTmpdir.never
    """keep temporary directory, {always, never, onerror}"""

    job_mult_memory: float = 1
    """memory multiplier"""

    job_mult_time: float = 1
    """running time multiplier"""

    merge_mult_memory: float = 1
    """memory multiplier for merging"""

    merge_mult_time: float = 1
    """running time multiplier for merging"""


class GatkPhaseByTransmission:
    de_novo_prior: float = 1e-8
    """use 1e-6 when interested in phasing de novos"""


class VariantPhasing(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    path_variant_annotation: Annotated[str, Field(examples=["../variant_annotation"])]

    tools_ngs_mapping: list[str] = []
    """expected tools for ngs mapping"""

    tools_variant_calling: list[str] = []
    """expected tools for variant calling"""

    phasings: list[str] = ["gatk_phasing_both"]

    ignore_chroms: list[str] = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*"]
    """patterns of chromosome names to ignore"""

    gatk_read_backed_phasing: GatkReadBackedPhasing = GatkReadBackedPhasing()

    gatk_phase_by_transmission: GatkPhaseByTransmission = GatkPhaseByTransmission()
