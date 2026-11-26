from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field

from snappy_pipeline.models import (
    RequiredBy,
    SizeString,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ToolModel,
)
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Exons,
    ExpectedPathSchema,
    ExpectedAlignments,  # noqa: F401  re-exported for the consumers of the alignments
    Features,
    Reference,
    TOOL,
)
from snappy_pipeline.workflows.common.strandedness import Strandedness
from snappy_pipeline.workflows.reference_index.model import ExpectedIndex


class ExpectedAlignmentQc(BaseModel):
    """Consumer-driven contract: the BAM QC reports of an ngs_mapping task that the VarFish
    export reads."""

    bamstats: str
    idxstats: str


class ExpectedRnaAlignments(BaseModel):
    """Consumer-driven contract: RNA alignments with their strandedness decision."""

    bam: str
    bai: str
    strandedness: str


class ExpectedGeneCounts(BaseModel):
    """Consumer-driven contract: the gene counts of an RNA library (STAR ``--quantMode``)."""

    gene_counts: str


class Tool(StrEnum):
    bwa = "bwa"
    bwa_mem2 = "bwa_mem2"
    bwa_mem3 = "bwa_mem3"
    minimap2 = "minimap2"
    star = "star"

    def is_dna(self):
        return self in {self.bwa, self.bwa_mem2, self.bwa_mem3, self.minimap2}


class NgsMappingDependsOn(SnappyModel):
    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
    ]
    """FASTQ source: an ``external_data`` task or an ``adapter_trimming`` task"""

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

    features: Annotated[Features, RequiredBy(Tool.star)] = ""

    exons: Annotated[Exons, RequiredBy(Tool.star)] = ""
    """Task with a BED file of ubiquitously expressed genes with a well-defined strand, from which
    STAR infers the strandedness"""


class BamCollectDoc(ToggleModel, ToolModel):
    window_length: Annotated[int, Field(gt=0)] = 1000


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
    """Additional arguments of ``bwa mem``"""


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


class NgsMapping(SnappyStepModel):
    bam_intermediate = True

    depends_on: NgsMappingDependsOn

    tool: Tool
    """Aligner to use for the NGS library"""

    bam_collect_doc: BamCollectDoc = BamCollectDoc()
    """Depth of coverage collection, mainly useful for genomes."""

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
