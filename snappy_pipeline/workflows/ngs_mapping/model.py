import enum
import os
from enum import Enum
from typing import Annotated

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from pydantic import DirectoryPath, Field, FilePath, field_validator, model_validator

from snappy_pipeline.models import EnumField, SizeString, SnappyModel, SnappyStepModel


class DnaMapper(Enum):
    BWA = "bwa"
    BWA2 = "bwa_mem2"


class LongDnaMapper(Enum):
    MINIMAP2 = "minimap2"


class RnaMapper(Enum):
    STAR = "star"


class Tools(SnappyModel):
    dna: Annotated[list[DnaMapper], EnumField(DnaMapper, [])]
    """Required if DNA analysis; otherwise, leave empty."""

    rna: Annotated[list[RnaMapper], EnumField(RnaMapper, [])]
    """Required if RNA analysis; otherwise, leave empty."""

    dna_long: Annotated[list[LongDnaMapper], EnumField(LongDnaMapper, [])]
    """Required if long-read mapper used; otherwise, leave empty."""


class TargetCoverageReportEntry(SnappyModel):
    """
    Mapping from enrichment kit to target region BED file, for either computing per--target
    region coverage or selecting targeted exons.

    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[os.PathLike, Field(examples=["path/to/targets.bed"])]


class TargetCoverageReport(SnappyModel):
    path_target_interval_list_mapping: list[TargetCoverageReportEntry] = []


class BamCollectDoc(SnappyModel):
    enabled: bool = False
    window_length: Annotated[int, Field(gt=0)] = 1000


class NgsChewFingerprint(SnappyModel):
    enabled: bool = True


class Bwa(SnappyModel):
    path_index: FilePath | str
    """Required if listed in ngs_mapping.tools.dna; otherwise, can be removed."""
    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    trim_adapters: bool = False
    mask_duplicates: bool = True
    split_as_secondary: bool = False
    """-M flag"""

    extra_flags: list[str] = []
    """ [ "-C" ] when molecular barcodes are processed with AGeNT in the somatic mode """

    @field_validator("path_index")
    def validate_bwa_path_index(cls, v):
        extensions = {".amb", ".ann", ".bwt", ".pac", ".sa"}
        prefix, ext = os.path.splitext(v)
        if ext:
            assert ext in extensions, f"unknown extension '{v}'"
        for extension in extensions:
            assert os.path.exists(prefix + extension), f"{v} does not exist"
        return prefix


class BwaMode(Enum):
    AUTO = "auto"
    BWA_ALN = "bwa-aln"
    BWA_MEM = "bwa-mem"


class BwaMem2(SnappyModel):
    path_index: FilePath | str
    """Required if listed in ngs_mapping.tools.dna; otherwise, can be removed."""

    bwa_mode: BwaMode = BwaMode.AUTO
    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    trim_adapters: bool = False
    mask_duplicates: bool = True
    split_as_secondary: bool = True
    """-M flag"""

    extra_flags: list[str] = []
    """[ "-C" ] when molecular barcodes are processed with AGeNT in the somatic mode"""

    @field_validator("path_index")
    def validate_bwa_mem2_path_index(cls, v):
        extensions = {".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"}
        prefix, ext = os.path.splitext(v)
        if ext:
            assert ext in extensions, f"unknown extension '{v}'"
        for extension in extensions:
            assert os.path.exists(prefix + extension), f"{v} does not exist"
        return prefix


class BarcodeTool(Enum):
    AGENT = "agent"


class Somatic(SnappyModel):
    mapping_tool: DnaMapper
    """Either bwa of bwa_mem2. The indices & other parameters are taken from mapper config"""

    barcode_tool: BarcodeTool = BarcodeTool.AGENT
    """Only agent currently implemented"""

    use_barcodes: bool = False
    recalibrate: bool = True


class Bqsr(SnappyModel):
    common_variants: str
    """Common germline variants (see /fast/work/groups/cubi/projects/biotools/static_data/app_support/GATK)"""


class AgentLibPrepType(Enum):
    HALO_PLEX = "halo"
    HALO_PLEX_HS = "hs"
    SURE_SELECT = "xt"
    SURE_SELECT_HS2 = "v2"
    SURE_SELECT_QXT = "qxt"


class AgentPrepare(SnappyModel):
    path: FilePath

    lib_prep_type: AgentLibPrepType = None
    """One of "halo" (HaloPlex), "hs" (HaloPlexHS), "xt" (SureSelect XT, XT2, XT HS), "v2" (SureSelect XT HS2) & "qxt" (SureSelect QXT)"""

    extra_args: list[str] = []
    """Consider "-polyG 8" for NovaSeq data & "-minFractionRead 50" for 100 cycles data"""


class AgentMarkDuplicatesConsensusMode(Enum):
    SINGLE = "SINGLE"
    HYBRID = "HYBRID"
    DUPLEX = "DUPLEX"


class AgentMarkDuplicates(SnappyModel):
    path: FilePath
    path_baits: FilePath
    consensus_mode: AgentMarkDuplicatesConsensusMode = None
    """One of "SINGLE", "HYBRID", "DUPLEX" """

    input_filter_args: list[str] = []
    """Consider -mm 13 (min base qual) -mr 13 (min barcode base qual) -mq 30 (min map qual)"""

    consensus_filter_args: list[str] = []

    extra_args: list[str] = []
    """Consider -d 1 (max nb barcode mismatch)"""


class Agent(SnappyModel):
    prepare: AgentPrepare
    mark_duplicates: AgentMarkDuplicates


class Star(SnappyModel):
    path_index: FilePath | str
    """Required if listed in ngs_mapping.tools.rna; otherwise, can be removed."""
    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    genome_load: str = "NoSharedMemory"
    raw_star_options: str = ""
    align_intron_max: int = 1000000  # ENCODE option
    align_intron_min: int = 20  # ENCODE option
    align_mates_gap_max: int = 1000000  # ENCODE option
    align_sjdb_overhang_min: int = 1  # ENCODE option
    align_sj_overhang_min: int = 8  # ENCODE option
    out_filter_mismatch_n_max: int = 999  # ENCODE option
    out_filter_mismatch_n_over_l_max: float = 0.04  # ENCODE option
    out_filter_multimap_n_max: int = 20  # ENCODE option
    out_filter_type: str = "BySJout"  # ENCODE option
    out_filter_intron_motifs: str = None
    """or for cufflinks: RemoveNoncanonical"""

    out_sam_strand_field: str = None
    """or for cufflinks: intronMotif"""

    transcriptome: bool = False
    """true to output transcript coordinate bam for RSEM"""

    trim_adapters: bool = False
    mask_duplicates: bool = False
    include_unmapped: bool = True

    @model_validator(mode="after")
    def ensure_star_index_files_exist(self):
        full_path = self.path_index
        # a lot of files should be in this dir, justtest these
        for indfile in ("Genome", "SA", "SAindex"):
            expected_path = os.path.join(full_path, indfile)
            if not os.path.exists(expected_path):  # pragma: no cover
                raise ValueError(f"Expected STAR index file {expected_path} does not exist!")
        return self


class Strand(enum.IntEnum):
    UNKNOWN = -1
    INFER = 0
    UNSTRANDED = 0
    FORWARD = 1
    REVERSE = 2


class Strandedness(SnappyModel):
    path_exon_bed: FilePath
    """Location of usually highly expressed genes. Known protein coding genes is a good choice"""

    strand: Strand = Strand.UNKNOWN
    """-1: unknown value, use infer_, 0: unstranded, 1: forward, 2: reverse (from featurecounts)"""

    threshold: float = 0.85
    """Minimum proportion of reads mapped to forward/reverse direction to call the protocol"""


class Minimap2(SnappyModel):
    mapping_threads: int = 16


class NgsMapping(SnappyStepModel):
    tools: Tools
    """Aligners to use for the different NGS library types"""

    path_link_in: DirectoryPath | str = ""
    """OPTIONAL Override data set configuration search paths for FASTQ files"""

    target_coverage_report: TargetCoverageReport | None = None
    """Thresholds for targeted sequencing coverage QC."""

    bam_collect_doc: BamCollectDoc = BamCollectDoc()
    """Depth of coverage collection, mainly useful for genomes."""

    ngs_chew_fingerprint: NgsChewFingerprint = NgsChewFingerprint()
    """Compute fingerprints with ngs-chew"""

    bwa: Bwa | None = None
    """Configuration for BWA"""

    bwa_mem2: BwaMem2 | None = None
    """Configuration for BWA-MEM2"""

    somatic: Somatic | None = None
    """
    Configuration for somatic ngs_calling
    (separate read groups, molecular barcodes & base quality recalibration)
    """

    bqsr: Bqsr | None = None

    agent: Agent | None = None

    star: Star | None = None
    """Configuration for STAR"""

    strandedness: Strandedness | None = None

    minimap2: Minimap2 | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for data_type in ("dna", "rna", "dna_long"):
            tool_list = getattr(self.tools, data_type)
            for tool in tool_list:
                if not getattr(self, tool):
                    raise ValueError(f"Tool {tool} not configured")
        return self
