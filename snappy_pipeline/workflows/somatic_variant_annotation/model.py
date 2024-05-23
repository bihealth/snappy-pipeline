import enum
from typing import Annotated

from pydantic import Field, model_validator

from models import SnappyStepModel, EnumField, SnappyModel


class Tool(enum.Enum):
    jannovar = "jannovar"
    vep = "vep"


class Dbnsfp(SnappyModel):
    col_contig: int = 1
    col_pos: int = 2
    columns: list[str] = []


class Jannovar(SnappyModel):
    path_jannovar_ser: str

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


class VepTxFlag(enum.Enum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class Vep(SnappyModel):
    cache_dir: str = ""
    """Defaults to $HOME/.vep Not a good idea on the cluster"""

    species: str = "homo_sapiens"

    assembly: str = "GRCh38"

    cache_version: int = 102
    """WARNING- this must match the wrapper's vep version!"""

    tx_flag: VepTxFlag = VepTxFlag.gencode_basic
    """The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged"."""

    pick_order: list[str] = [
        "biotype",
        "mane",
        "appris",
        "tsl",
        "ccds",
        "canonical",
        "rank",
        "length",
    ]
    num_threads: int = 8
    buffer_size: int = 1000
    output_options: list[str] = ["everything"]


class SomaticVariantAnnotation(SnappyStepModel):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.jannovar, Tool.vep])]

    path_somatic_variant_calling: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    tools_ngs_mapping: list[str] = []
    """default to those configured for ngs_mapping"""

    tools_somatic_variant_calling: list[str] = []
    """default to those configured for somatic_variant_calling"""

    jannovar: Jannovar | None = None

    vep: Vep | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
