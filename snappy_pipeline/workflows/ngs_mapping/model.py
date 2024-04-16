import re
import typing
from enum import Enum
from typing import Annotated, Self

from annotated_types import Predicate
from pydantic import BaseModel, conint, ConfigDict, model_validator, model_serializer

size_string_regexp = re.compile(r"[. 0-9]+([KMGTP])")
SizeString = Annotated[str, Predicate(lambda s: size_string_regexp.match(s) is not None)]


class DnaMapper(Enum):
    BWA = "bwa"
    BWA2 = "bwa_mem2"


class LongDnaMapper(Enum):
    MINIMAP2 = "minimap2"


class RnaMapper(Enum):
    STAR = "star"


class Tools(BaseModel):
    dna: list[DnaMapper] = []
    """Required if DNA analysis; otherwise, leave empty. Example: 'bwa'."""

    rna: list[RnaMapper] = []
    """Required if RNA analysis; otherwise, leave empty. Example: 'star'."""

    dna_long: list[LongDnaMapper] = []
    """Required if long-read mapper used; otherwise, leave empty. Example: 'minimap2'."""


class TargetCoverageReportEntry(BaseModel):
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

    name: str
    """'IDT_xGen_V1_0'"""
    pattern: str
    """'xGen Exome Research Panel V1\\.0*'"""
    path: str
    """'path/to/targets.bed'"""


class TargetCoverageReport(BaseModel):
    path_target_interval_list_mapping: list[TargetCoverageReportEntry] = []


class BamCollectDoc(BaseModel):
    enabled: bool = False
    window_length: conint(ge=0) = 1000


class NgsChewFingerprint(BaseModel):
    enabled: bool = True


class Bwa(BaseModel):
    path_index: str = None
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


class BwaMode(Enum):
    AUTO = "auto"
    BWA_ALN = "bwa-aln"
    BWA_MEM = "bwa-mem"


class BwaMem2(BaseModel):
    path_index: str = None
    """Required if listed in ngs_mapping.tools.dna; otherwise, can be removed."""

    bwa_mode: BwaMode = BwaMode.AUTO
    num_threads_align: int | None = 16
    num_threads_trimming: int | None = 8
    num_threads_bam_view: int | None = 4
    num_threads_bam_sort: int | None = 4
    memory_bam_sort: SizeString | None = "4G"
    trim_adapters: bool | None = False
    mask_duplicates: bool | None = True
    split_as_secondary: bool | None = True
    """-M flag"""

    extra_flags: list[str] = []
    """[ "-C" ] when molecular barcodes are processed with AGeNT in the somatic mode"""


class BarcodeTool(Enum):
    AGENT = "agent"


class Somatic(BaseModel):
    mapping_tool: DnaMapper = None
    """Either bwa of bwa_mem2. The indices & other parameters are taken from mapper config"""

    barcode_tool: BarcodeTool = BarcodeTool.AGENT
    """Only agent currently implemented"""

    use_barcodes: bool = False
    recalibrate: bool = True


class Bqsr(BaseModel):
    common_variants: str
    """Common germline variants (see /fast/work/groups/cubi/projects/biotools/static_data/app_support/GATK)"""


class AgentLibPrepType(Enum):
    HALO_PLEX = "halo"
    HALO_PLEX_HS = "hs"
    SURE_SELECT = "xt"
    SURE_SELECT_HS2 = "v2"
    SURE_SELECT_QXT = "qxt"


class AgentPrepare(BaseModel):
    path: str = None
    lib_prep_type: AgentLibPrepType = None
    """One of "halo" (HaloPlex), "hs" (HaloPlexHS), "xt" (SureSelect XT, XT2, XT HS), "v2" (SureSelect XT HS2) & "qxt" (SureSelect QXT)"""

    extra_args: list[str] = []
    """Consider "-polyG 8" for NovaSeq data & "-minFractionRead 50" for 100 cycles data"""


class AgentMarkDuplicatesConsensusMode(Enum):
    SINGLE = "SINGLE"
    HYBRID = "HYBRID"
    DUPLEX = "DUPLEX"


class AgentMarkDuplicates(BaseModel):
    path: str = None
    path_baits: str = None
    consensus_mode: AgentMarkDuplicatesConsensusMode = None
    """One of "SINGLE", "HYBRID", "DUPLEX" """

    input_filter_args: list[str] = []
    """Consider -mm 13 (min base qual) -mr 13 (min barcode base qual) -mq 30 (min map qual)"""

    consensus_filter_args: list[str] = []

    extra_args: list[str] = []
    """Consider -d 1 (max nb barcode mismatch)"""


class Agent(BaseModel):
    prepare: AgentPrepare
    mark_duplicates: AgentMarkDuplicates


class Star(BaseModel):
    path_index: str
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
    out_filter_intron_motifs: str | None = None
    """or for cufflinks: RemoveNoncanonical"""

    out_sam_strand_field: str | None = None
    """or for cufflinks: intronMotif"""

    transcriptome: bool = False
    """true to output transcript coordinate bam for RSEM"""

    trim_adapters: bool = False
    mask_duplicates: bool = False
    include_unmapped: bool = True


class Strand(Enum):
    UNKNOWN = -1
    INFER = 0
    UNSTRANDED = 0
    FORWARD = 1
    REVERSE = 2


class Strandedness(BaseModel):
    path_exon_bed: str
    """Location of usually highly expressed genes. Known protein coding genes is a good choice"""

    strand: Strand = Strand.UNKNOWN
    """-1: unknown value, use infer_, 0: unstranded, 1: forward, 2: reverse (from featurecounts)"""

    threshold: float = 0.85
    """Minimum proportion of reads mapped to forward/reverse direction to call the protocol"""


class Minimap2(BaseModel):
    mapping_threads: int = 16


class NgsMapping(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
    )

    tools: Tools = Tools()
    """Aligners to use for the different NGS library types"""

    path_link_in: str | None = None
    """OPTIONAL Override data set configuration search paths for FASTQ files"""

    target_coverage_report: TargetCoverageReport | None = TargetCoverageReport()
    """Thresholds for targeted sequencing coverage QC."""

    bam_collect_doc: BamCollectDoc | None = BamCollectDoc()
    """Depth of coverage collection, mainly useful for genomes."""

    ngs_chew_fingerprint: NgsChewFingerprint | None = NgsChewFingerprint()
    """Compute fingerprints with ngs-chew"""

    bwa: Bwa | None = Bwa()
    """Configuration for BWA"""

    bwa_mem2: BwaMem2 | None = BwaMem2()
    """Configuration for BWA-MEM2"""

    somatic: Somatic | None = Somatic()
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
                print(tool)
                if tool not in (tool_config := getattr(self, tool)):
                    raise ValueError(f"Tool {tool} not configured in {tool_config}")
        return self

    @model_serializer(mode="wrap")
    def _serialize(self, handler):
        d = handler(self)
        print(d)
        return d
