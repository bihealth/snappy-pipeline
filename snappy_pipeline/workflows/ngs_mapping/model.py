import enum
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import (
    SizeString,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ToolModel,
    validators,
)
from snappy_pipeline.workflows.abstract.protocol import (
    CommonVariants,
    DataSignature,
    DataType,
    Exons,
    ExpectedPathSchema,
    Features,
    Kits,
    Reference,
    TOOL,
)
from snappy_pipeline.workflows.reference_index.model import ExpectedIndex


class ExpectedAlignments(BaseModel):
    """Consumer-driven contract: expected output keys from an ngs_mapping upstream task."""

    bam: str
    bai: str


class ExpectedAlignmentQc(BaseModel):
    """Consumer-driven contract: the BAM QC reports of an ngs_mapping task."""

    bamstats: str
    flagstats: str
    idxstats: str
    alfred_qc: str | None = None
    """only if the task enables ``target_coverage_report``"""


class ExpectedRnaAlignments(BaseModel):
    """Consumer-driven contract: RNA alignments with their strandedness decision."""

    bam: str
    bai: str
    strandedness: str


class ExpectedGeneCounts(BaseModel):
    """Consumer-driven contract: the gene counts of an RNA library (STAR ``--quantMode``)."""

    gene_counts: str


class NgsMappingDependsOn(SnappyModel):
    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
    ]
    """FASTQ source: ``data_sets``, an ``external_data`` task or an ``adapter_trimming`` task"""

    index: Annotated[
        str,
        DataSignature(
            DataType.INDEX,
            frozenset({TOOL, ("dna", "rna")}),
        ),
        ExpectedPathSchema(ExpectedIndex),
    ]
    """Task with the index of the mapper, tagged with its tool (mbcs: its ``mapping_tool``)."""

    reference: Reference

    features: Features = ""

    kits: Kits = ""

    exons: Exons = ""
    """Task with a BED file of ubiquitously expressed genes with a well-defined strand; needed
    for star, which infers the strandedness from it."""

    common_variants: CommonVariants = ""
    """Needed by ``tool: mbcs`` with ``recalibrate`` (BQSR)"""


class DnaMapper(StrEnum):
    BWA = "bwa"
    BWA_MEM2 = "bwa_mem2"


class LongDnaMapper(StrEnum):
    MINIMAP2 = "minimap2"


class RnaMapper(StrEnum):
    STAR = "star"


class MetaTool(StrEnum):
    MBCS = "mbcs"


class Tool(StrEnum):
    bwa = "bwa"
    bwa_mem2 = "bwa_mem2"
    minimap2 = "minimap2"
    star = "star"
    mbcs = "mbcs"

    def is_dna(self):
        return self in {self.bwa, self.bwa_mem2, self.minimap2, self.mbcs}

    def is_rna(self):
        return self in {self.star}

    def supports_long_reads(self):
        return self in {self.bwa_mem2, self.minimap2}

    def get_tags(self):
        tags = set()
        if self.is_dna():
            tags |= {"dna"}
        if self.is_rna():
            tags |= {"rna"}
        if self.supports_long_reads():
            tags |= {"long_read"}
        if self == self.mbcs:
            tags |= {"mbcs", "meta"}
        return tags


class TargetCoverageReport(ToggleModel, ToolModel):
    """The report takes the ``targets`` of each library's kit from ``depends_on.kits``."""


class BamCollectDoc(ToggleModel, ToolModel):
    window_length: Annotated[int, Field(gt=0)] = 1000


class NgsChewFingerprint(ToggleModel, ToolModel):
    pass


class BwaMode(StrEnum):
    AUTO = "auto"
    BWA_ALN = "bwa-aln"
    BWA_MEM = "bwa-mem"


class BwaMapper(ToolModel):
    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    trim_adapters: bool = False
    mask_duplicates: bool = True

    split_as_secondary: bool = False
    """-M flag"""

    platform: str = "ILLUMINA"
    """Platform (``PL``) of the read group"""

    extra_arguments: list[str] = []
    """Additional arguments of ``bwa mem``, such as ``["-C"]`` when AGeNT processes molecular
    barcodes (``mbcs``)"""


class Bwa(BwaMapper):
    aln_max_read_length: int = 75
    """Libraries with an average read length up to this value (estimated from the first 100000
    reads) are mapped with ``bwa aln``, longer ones with ``bwa mem``"""


class BwaMem2(BwaMapper):
    pass


class BarcodeTool(StrEnum):
    AGENT = "agent"


class AgentLibPrepType(StrEnum):
    HALO_PLEX = "halo"
    HALO_PLEX_HS = "hs"
    SURE_SELECT = "xt"
    SURE_SELECT_HS2 = "v2"
    SURE_SELECT_QXT = "qxt"


class AgentPrepare(ToolModel):
    path: ConfigPath

    lib_prep_type: AgentLibPrepType = None
    """One of "halo" (HaloPlex), "hs" (HaloPlexHS), "xt" (SureSelect XT, XT2, XT HS), "v2" (SureSelect XT HS2) & "qxt" (SureSelect QXT)"""

    extra_arguments: list[str] = []
    """Consider "-polyG 8" for NovaSeq data & "-minFractionRead 50" for 100 cycles data"""


class AgentMarkDuplicatesConsensusMode(enum.StrEnum):
    SINGLE = "SINGLE"
    HYBRID = "HYBRID"
    DUPLEX = "DUPLEX"


class AgentMarkDuplicates(ToolModel):
    path: ConfigPath
    consensus_mode: AgentMarkDuplicatesConsensusMode = None
    """One of "SINGLE", "HYBRID", "DUPLEX" """

    input_filter_args: list[str] = []
    """Consider -mm 13 (min base qual) -mr 13 (min barcode base qual) -mq 30 (min map qual)"""

    consensus_filter_args: list[str] = []

    extra_arguments: list[str] = []
    """Consider -d 1 (max nb barcode mismatch)"""


class Agent(SnappyModel):
    prepare: AgentPrepare
    mark_duplicates: AgentMarkDuplicates


class Star(ToolModel):
    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    genome_load: str = "NoSharedMemory"
    raw_star_options: str = ""
    platform: str = "ILLUMINA"
    """Platform (``PL``) of the read group"""

    align_intron_max: int = 1000000  # ENCODE option
    align_intron_min: int = 20  # ENCODE option
    align_mates_gap_max: int = 1000000  # ENCODE option
    align_sjdb_overhang_min: int = 1  # ENCODE option
    align_sj_overhang_min: int = 8  # ENCODE option
    out_filter_mismatch_n_max: int = 999  # ENCODE option
    out_filter_mismatch_n_over_l_max: float = 0.04  # ENCODE option
    out_filter_multimap_n_max: int = 20  # ENCODE option
    out_filter_type: str = "BySJout"  # ENCODE option
    out_filter_intron_motifs: str = "None"
    """or for cufflinks: RemoveNoncanonical"""

    out_sam_strand_field: str = "None"
    """or for cufflinks: intronMotif"""

    transcriptome: bool = False
    """true to output transcript coordinate bam for RSEM"""

    trim_adapters: bool = False
    mask_duplicates: bool = False
    include_unmapped: bool = True


class Strand(enum.IntEnum):
    UNKNOWN = -1
    INFER = 0
    UNSTRANDED = 0
    FORWARD = 1
    REVERSE = 2


class Strandedness(ToolModel):
    strand: Strand = Strand.UNKNOWN
    """-1: unknown value, use infer_, 0: unstranded, 1: forward, 2: reverse (from featurecounts)"""

    threshold: float = 0.85
    """Minimum proportion of reads mapped to forward/reverse direction to call the protocol"""


class Minimap2(ToolModel):
    mapping_threads: int = 16


class Mbcs(SnappyModel):
    mapping_tool: DnaMapper
    """Either bwa or bwa_mem2. The indices & other parameters are taken from mapper config"""

    barcode_tool: BarcodeTool = BarcodeTool.AGENT
    """Only agent currently implemented"""

    use_barcodes: bool = False
    recalibrate: bool = True


TOOL_DEPENDENCIES = {Tool.star: ("features", "exons")}


class NgsMapping(SnappyStepModel):
    depends_on: NgsMappingDependsOn

    tool: Tool
    """Aligner to use for the NGS library"""

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

    agent: Agent | None = None

    star: Star | None = None
    """Configuration for STAR"""

    strandedness: Strandedness = Strandedness()

    minimap2: Minimap2 | None = None

    mbcs: Mbcs | None = None
    """
    Configuration for somatic ngs_calling
    (separate read groups, molecular barcodes & base quality recalibration)
    """

    def tool_tag(self) -> str | None:
        """mbcs maps with ``mbcs.mapping_tool``, whose index it reads."""
        if self.tool == Tool.mbcs and self.mbcs:
            return self.mbcs.mapping_tool
        return self.tool

    @model_validator(mode="after")
    def check_mbcs_prerequisites(self):
        if self.mbcs:
            tool = self.mbcs.mapping_tool
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")

            if self.mbcs.use_barcodes:
                if not self.agent:
                    raise ValueError("Agent configuration required for MBCS")
                if not self.depends_on.kits:
                    raise ValueError("mbcs with use_barcodes needs depends_on.kits")

            if self.mbcs.recalibrate:
                if not self.depends_on.common_variants:
                    raise ValueError("mbcs with recalibrate needs depends_on.common_variants")
        return self

    @model_validator(mode="after")
    def check_target_coverage_report_kits(self):
        if self.target_coverage_report and self.target_coverage_report.enabled:
            if not self.depends_on.kits:
                raise ValueError("target_coverage_report needs depends_on.kits")
        return self

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
