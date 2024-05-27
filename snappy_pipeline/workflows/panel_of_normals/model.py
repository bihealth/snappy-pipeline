import enum
from typing import Annotated, Self

from pydantic import Field, model_validator, DirectoryPath

from snappy_pipeline.models import SnappyStepModel, EnumField, SnappyModel, KeepTmpdir


class Tool(enum.StrEnum):
    mutect2 = "mutect2"


class Mutect2(SnappyModel):
    path_normals_list: str = ""

    germline_resource: str

    java_options: str = " -Xmx16g"

    num_cores: int = 2
    """number of cores to use locally"""

    window_length: int = 100000000
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


class CnvKit(SnappyModel):
    path_normals_list: str = ""
    """Optional file listing libraries to include in panel"""

    path_target_regions: str = ""
    """Bed files of targetted regions (Missing when creating a panel of normals for WGS data)"""

    access: str = ""
    """Access bed file (output/cnvkit.access/out/cnvkit.access.bed when create_cnvkit_acces was run)"""

    annotate: str = ""
    """[target] Optional targets annotations"""

    target_avg_size: int = 0
    """[target] Average size of split target bins (0: use default value)"""

    bp_per_bin: int = 50000
    """[autobin] Expected base per bin"""

    split: bool = True
    """[target] Split large intervals into smaller ones"""

    antitarget_avg_size: int = 0
    """[antitarget] Average size of antitarget bins (0: use default value)"""

    min_size: int = 0
    """[antitarget] Min size of antitarget bins (0: use default value)"""

    min_mapq: int = 0
    """[coverage] Mininum mapping quality score to count a read for coverage depth"""

    count: bool = False
    """[coverage] Alternative couting algorithm"""

    min_cluster_size: int = 0
    """[reference] Minimum cluster size to keep in reference profiles. 0 for no clustering"""

    gender: str = ""
    """[reference] Specify the chromosomal sex of all given samples as male or female. Guess when missing"""

    male_reference: bool = False
    """[reference & sex] Create male reference"""

    gc_correction: bool = True
    """[reference] Use GC correction"""

    edge_correction: bool = True
    """[reference] Use edge correction"""

    rmask_correction: bool = True
    """[reference] Use rmask correction"""

    drop_low_coverage: bool = False
    """[metrics] Drop very-low-coverage bins before calculations"""


class Access(SnappyModel):
    """Creates access file for cnvkit, based on genomic sequence & excluded regions (optionally)"""

    exclude: list[str] = []
    """[access] Bed file of regions to exclude (mappability, blacklisted, ...)"""

    min_gap_size: int = 0
    """[access] Minimum gap size between accessible sequence regions (0: use default value)"""


class GenomeName(enum.StrEnum):
    hg18 = "hg18"
    hg19 = "hg19"
    hg38 = "hg38"
    mm9 = "mm9"
    mm10 = "mm10"
    rn4 = "rn4"
    rn5 = "rn5"
    rn6 = "rn6"
    canFam3 = "canFam3"


class PureCn(SnappyModel):
    path_normals_list: str = ""
    """Optional file listing libraries to include in panel"""

    path_bait_regions: str
    """
    Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect),
    recommended by PureCN author
    """

    path_genomicsDB: str
    """Mutect2 genomicsDB created during panel_of_normals"""

    genome_name: Annotated[GenomeName, EnumField(GenomeName)]

    enrichment_kit_name: str = "unknown"
    """For filename only..."""

    mappability: str = ""
    """
    GRCh38:
    /fast/work/groups/cubi/projects/biotools/static_data/app_support/PureCN/hg38/mappability.bw
    """

    reptiming: str = ""
    """Nothing for GRCh38"""

    seed: int = 1234567


class PanelOfNormals(SnappyStepModel):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mutect2], min_length=1)]

    path_ngs_mapping: DirectoryPath | str

    ignore_chroms: Annotated[
        list[str],
        Field(
            examples=[
                "NC_007605",
                "hs37d5",
                "chrEBV",
                "*_decoy",
                "HLA-*",
                "GL000220.*",
                "chrEBV",
                "HPV*",
                "CMV",
                "HBV",
                "HCV-*",
                "HIV-*",
                "KSHV",
                "HTLV-1",
                "MCV",
                "*_decoy",
                "chrUn_GL00220*",
                "SV40",
            ]
        ),
    ] = []
    """Patterns of contig names to ignore"""

    mutect2: Mutect2 | None = None

    cnvkit: CnvKit | None = None

    access: Access = Access()

    purecn: PureCn | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
