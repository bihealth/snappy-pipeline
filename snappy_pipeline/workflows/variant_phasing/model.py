from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import KeepTmpdir, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedAnnotatedGermlineVariants


class GatkReadBackedPhasing(SnappyModel):
    phase_quality_threshold: float = 20.0
    """quality threshold for phasing"""

    window_length: int = 5000000
    """split input into windows of this size, each triggers a job"""

    num_jobs: int = 1000
    """number of windows to process in parallel"""

    use_profile: bool = True
    """use Snakemake profile for parallel processing"""

    restart_runtimes: int = 0
    """number of runtimes to re-launch jobs in case of failure"""

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

    job_mult_runtime: float = 1
    """running runtime multiplier"""

    merge_mult_memory: float = 1
    """memory multiplier for merging"""

    merge_mult_runtime: float = 1
    """running runtime multiplier for merging"""


class GatkPhaseByTransmission(SnappyModel):
    de_novo_prior: float = 1e-8
    """use 1e-6 when interested in phasing de novos"""


class ExpectedPhasedVariants(SnappyModel):
    """Consumer-driven contract: expected output keys from variant_phasing."""

    vcf: str
    vcf_tbi: str


class VariantPhasingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    variant_annotation: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "annotated"})),
        ExpectedPathSchema(ExpectedAnnotatedGermlineVariants),
    ] = "variant_annotation"


class VariantPhasing(SnappyStepModel):
    depends_on: VariantPhasingDependsOn = Field(default_factory=VariantPhasingDependsOn)

    phasings: list[str] = ["gatk_phasing_both"]

    ignore_chroms: list[str] = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*"]
    """patterns of chromosome names to ignore"""

    gatk_read_backed_phasing: GatkReadBackedPhasing = GatkReadBackedPhasing()

    gatk_phase_by_transmission: GatkPhaseByTransmission = GatkPhaseByTransmission()
