import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.annotation import Vep


class Tool(enum.StrEnum):
    jannovar = "jannovar"
    vep = "vep"


class Dbnsfp(SnappyModel):
    col_contig: int = 1
    col_pos: int = 2
    columns: list[str] = []


class Jannovar(SnappyModel):
    path_jannovar_ser: str
    """Path to serialized Jannovar database"""

    flag_off_target: bool

    dbnsfp: Dbnsfp
    """configuration for default genome release, needs change if differing"""

    annotation_tracks_bed: list[str] = []

    annotation_tracks_tsv: list[str] = []

    annotation_tracks_vcf: list[str] = []

    window_length: int = 50000000
    """split input into windows of this size, each triggers a job"""

    num_jobs: int = 100
    """number of windows to process in parallel"""

    use_profile: bool = True
    """use Snakemake profile for parallel processing"""

    restart_times: int = 5
    """number of times to re-launch jobs in case of failure"""

    max_jobs_per_second: int = 10
    """throttling of job creation"""

    max_status_checks_per_second: int = 10
    """throttling of status checks"""

    ignore_chroms: list[str] = ["NC_007605", "hs37d5", "chrEBV", "GL*", "*_decoy", "HLA-*"]
    """patterns of chromosome names to ignore"""


class SomaticVariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.jannovar, Tool.vep], min_length=1)]

    path_somatic_variant_calling: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    jannovar: Jannovar | None = None

    vep: Vep | None = None
