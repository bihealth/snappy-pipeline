import enum
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import (
    SizeString,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ToolModel,
    validators,
)
from snappy_pipeline.workflows.abstract.protocol import (
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
    """Task with the index of the mapper, tagged with its tool."""

    reference: Reference

    features: Features = ""

    kits: Kits = ""

    exons: Exons = ""
    """Task with a BED file of ubiquitously expressed genes with a well-defined strand; needed
    for star, which infers the strandedness from it."""


class Tool(StrEnum):
    bwa = "bwa"
    bwa_mem2 = "bwa_mem2"
    bwa_mem3 = "bwa_mem3"
    minimap2 = "minimap2"
    star = "star"

    def is_dna(self):
        return self in {self.bwa, self.bwa_mem2, self.bwa_mem3, self.minimap2}

    def is_rna(self):
        return self in {self.star}

    def supports_long_reads(self):
        return self in {self.bwa_mem2, self.bwa_mem3, self.minimap2}

    def get_tags(self):
        tags = set()
        if self.is_dna():
            tags |= {"dna"}
        if self.is_rna():
            tags |= {"rna"}
        if self.supports_long_reads():
            tags |= {"long_read"}
        return tags


class TargetCoverageReport(ToggleModel, ToolModel):
    """The report takes the ``targets`` of each library's kit from ``depends_on.kits``."""

    ignore_read_groups: bool = True
    """Compute one set of metrics per library, not per read group (``alfred qc --ignore``)"""


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

    bases_per_batch: Annotated[int, Field(gt=0)] = 100_000_000
    """Input bases per batch of ``bwa mem`` (``-K``). A fixed batch makes the alignments
    independent of ``num_threads_align``: without it, the batch is 10 Mbp per thread, and the
    insert size distribution is estimated per batch."""

    platform: str = "ILLUMINA"
    """Platform (``PL``) of the read group"""

    extra_arguments: list[str] = []
    """Additional arguments of ``bwa mem``, such as ``["-C"]`` for reads that AGeNT trimmed, so that
    their molecular barcodes become tags for AGeNT's duplicate marking"""


class Bwa(BwaMapper):
    aln_max_read_length: int = 75
    """Libraries with an average read length up to this value (estimated from the first 100000
    reads) are mapped with ``bwa aln``, longer ones with ``bwa mem``"""


class BwaMem2(BwaMapper):
    pass


class BwaMem3(BwaMapper):
    """BWA-MEM3 (Fulcrum Genomics) with its default settings: the alignments of BWA-MEM2, plus
    the tags ``MQ`` (mate mapping quality) and ``HN`` (number of hits)"""


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
    platform: str = "PACBIO"
    """Platform (``PL``) of the read group"""

    soft_clip_supplementary: bool = True
    """Use soft clipping for supplementary alignments (``-Y``)"""

    output_md_tag: bool = True
    """Write the ``MD`` tag (``--MD``)"""

    memory_bam_sort: SizeString = "4G"
    """Memory per thread of ``samtools sort`` (``-m``)"""


TOOL_DEPENDENCIES = {Tool.star: ("features", "exons")}


class NgsMapping(SnappyStepModel):
    bam_intermediate = True

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

    bwa_mem3: BwaMem3 | None = None
    """Configuration for BWA-MEM3"""

    star: Star | None = None
    """Configuration for STAR"""

    strandedness: Strandedness = Strandedness()

    minimap2: Minimap2 | None = None

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
