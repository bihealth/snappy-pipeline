import enum
import os
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import (
    SizeString,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ResolvablePathPrefix,
    ResolvablePath,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.adapter_trimming.model import ExpectedTrimmedRawFastq
from snappy_pipeline.workflows.link_in.model import ExpectedLinkedRawFastq
from snappy_pipeline.workflows.reference_index.model import ExpectedReferenceIndexFiles


class ExpectedAlignments(BaseModel):
    """Consumer-driven contract: expected output keys from an ngs_mapping upstream task."""

    bam: str
    bai: str


class NgsMappingDependsOn(SnappyModel):
    # External FASTQ source. Usually points to a dedicated link_in task.
    link_in: Annotated[
        str,
        DataSignature(DataType.RAW),
        ExpectedPathSchema(ExpectedLinkedRawFastq),
    ] = ""

    # Optional in-pipeline FASTQ source, e.g. adapter_trimming output.
    adapter_trimming: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"trimmed"})),
        ExpectedPathSchema(ExpectedTrimmedRawFastq),
    ] = ""

    # Optional upstream index provider task.
    reference_index: Annotated[
        str,
        DataSignature(
            DataType.INDEX,
            frozenset({("bwa", "bwa_mem2", "minimap2", "star"), ("dna", "rna")}),
        ),
        ExpectedPathSchema(ExpectedReferenceIndexFiles),
    ] = ""


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

    path: ResolvablePath = Field(examples=["path/to/targets.bed"])


class TargetCoverageReport(ToggleModel):
    path_target_interval_list_mapping: list[TargetCoverageReportEntry] = []


class BamCollectDoc(ToggleModel):
    window_length: Annotated[int, Field(gt=0)] = 1000


class NgsChewFingerprint(ToggleModel):
    pass


class BwaMode(StrEnum):
    AUTO = "auto"
    BWA_ALN = "bwa-aln"
    BWA_MEM = "bwa-mem"


class BwaMapper(SnappyModel):
    path_index: ResolvablePathPrefix
    """Path prefix for BWA index files (e.g., "path/to/GRCh38" without ".amb" extension)"""

    num_threads_align: int = 16
    num_threads_trimming: int = 8
    num_threads_bam_view: int = 4
    num_threads_bam_sort: int = 4
    memory_bam_sort: SizeString = "4G"
    trim_adapters: bool = False
    mask_duplicates: bool = True

    split_as_secondary: bool = False
    """-M flag"""

    extra_args: list[str] = []
    """[ "-C" ] when molecular barcodes are processed with AGeNT in the somatic mode"""


class Bwa(BwaMapper):
    @model_validator(mode="after")
    def validate_bwa_path_index(self):
        import logging

        v = self.path_index
        extensions = {".amb", ".ann", ".bwt", ".pac", ".sa"}
        prefix, ext = os.path.splitext(v)
        if ext:
            if ext in {".fa", ".fasta"}:
                prefix += ext
            else:
                if ext not in extensions:
                    logging.warning(f"unknown extension '{v}'")
        for extension in extensions:
            sidecar = prefix + extension
            alt_sidecar = os.path.splitext(prefix)[0] + extension
            if not (os.path.exists(sidecar) or os.path.exists(alt_sidecar)):
                logging.warning(f"missing BWA index sidecar file: {sidecar} (or {alt_sidecar})")
        self.path_index = prefix
        return self


class BwaMem2(BwaMapper):
    @model_validator(mode="after")
    def validate_bwa_mem2_path_index(self):
        import logging

        v = self.path_index
        extensions = {".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"}
        prefix, ext = os.path.splitext(v)
        if ext:
            if ext in {".fa", ".fasta"}:
                prefix += ext
            else:
                if ext not in extensions:
                    logging.warning(f"unknown extension '{v}'")
        for extension in extensions:
            sidecar = prefix + extension
            alt_sidecar = os.path.splitext(prefix)[0] + extension
            if not (os.path.exists(sidecar) or os.path.exists(alt_sidecar)):
                logging.warning(
                    f"missing BWA-MEM2 index sidecar file: {sidecar} (or {alt_sidecar})"
                )
        self.path_index = prefix
        return self


class BarcodeTool(StrEnum):
    AGENT = "agent"


class Bqsr(SnappyModel):
    common_variants: ResolvablePath
    """Common germline variants (see /fast/work/groups/cubi/projects/biotools/static_data/app_support/GATK)"""


class AgentLibPrepType(StrEnum):
    HALO_PLEX = "halo"
    HALO_PLEX_HS = "hs"
    SURE_SELECT = "xt"
    SURE_SELECT_HS2 = "v2"
    SURE_SELECT_QXT = "qxt"


class AgentPrepare(SnappyModel):
    path: str

    lib_prep_type: AgentLibPrepType = None
    """One of "halo" (HaloPlex), "hs" (HaloPlexHS), "xt" (SureSelect XT, XT2, XT HS), "v2" (SureSelect XT HS2) & "qxt" (SureSelect QXT)"""

    extra_args: list[str] = []
    """Consider "-polyG 8" for NovaSeq data & "-minFractionRead 50" for 100 cycles data"""


class AgentMarkDuplicatesConsensusMode(enum.StrEnum):
    SINGLE = "SINGLE"
    HYBRID = "HYBRID"
    DUPLEX = "DUPLEX"


class AgentMarkDuplicates(SnappyModel):
    path: str
    path_baits: str
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
    path_index: str
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
    out_filter_intron_motifs: str = "None"
    """or for cufflinks: RemoveNoncanonical"""

    out_sam_strand_field: str = "None"
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
    path_exon_bed: str
    """Location of usually highly expressed genes. Known protein coding genes is a good choice"""

    strand: Strand = Strand.UNKNOWN
    """-1: unknown value, use infer_, 0: unstranded, 1: forward, 2: reverse (from featurecounts)"""

    threshold: float = 0.85
    """Minimum proportion of reads mapped to forward/reverse direction to call the protocol"""


class Minimap2(SnappyModel):
    mapping_threads: int = 16

    path_index: str


class Mbcs(SnappyModel):
    mapping_tool: DnaMapper
    """Either bwa or bwa_mem2. The indices & other parameters are taken from mapper config"""

    barcode_tool: BarcodeTool = BarcodeTool.AGENT
    """Only agent currently implemented"""

    use_barcodes: bool = False
    recalibrate: bool = True


class NgsMapping(SnappyStepModel):
    depends_on: NgsMappingDependsOn = Field(default_factory=NgsMappingDependsOn)

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

    bqsr: Bqsr | None = None

    agent: Agent | None = None

    star: Star | None = None
    """Configuration for STAR"""

    strandedness: Strandedness | None = None

    minimap2: Minimap2 | None = None

    mbcs: Mbcs | None = None
    """
    Configuration for somatic ngs_calling
    (separate read groups, molecular barcodes & base quality recalibration)
    """

    @model_validator(mode="after")
    def check_mbcs_prerequisites(self):
        if self.mbcs:
            tool = self.mbcs.mapping_tool
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")

            if self.mbcs.use_barcodes:
                if not self.agent:
                    raise ValueError("Agent configuration required for MBCS")

            if self.mbcs.recalibrate:
                if not self.bqsr:
                    raise ValueError("BQSR configuration required for MBCS")
        return self
