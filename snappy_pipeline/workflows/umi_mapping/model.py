from enum import StrEnum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
    Umis,
)
from snappy_pipeline.workflows.reference_index.model import ExpectedIndex


class Tool(StrEnum):
    fgumi = "fgumi"


class Consensus(StrEnum):
    duplex = "duplex"
    simplex = "simplex"


class UmiMappingDependsOn(SnappyModel):
    reads: Annotated[str, DataSignature(DataType.RAW)]
    """FASTQ source: an ``external_data`` task or an ``adapter_trimming`` task"""

    index: Annotated[
        str,
        DataSignature(DataType.INDEX, frozenset({"bwa_mem3", "dna"})),
        ExpectedPathSchema(ExpectedIndex),
    ]
    """Task with the bwa-mem3 index of the genome: a ``reference_index`` task with
    ``tool: bwa_mem3``, or an ``external_data`` task with
    ``produces: {type: index, tags: [bwa_mem3, dna]}``"""

    reference: Reference

    umis: Annotated[Umis, RequiredBy(Tool.fgumi, when="fgumi.correct")] = ""
    """Task with the fixed UMIs of the kit, one sequence per line, against which ``fgumi correct``
    corrects the extracted UMIs (with the ``correct`` section)"""


class Aligner(SnappyModel):
    """``bwa-mem3 mem``, which aligns the raw reads and realigns the consensus reads"""

    num_threads: Annotated[int, Field(gt=0)] = 16
    """Threads of bwa-mem3, which are the threads of the grouping job and of the consensus job"""

    bases_per_batch: Annotated[int, Field(gt=0)] = 150_000_000
    """Input bases per batch of the raw-read alignment (``-K``), which makes the alignments
    independent of the threads"""

    extra_arguments: list[str] = []
    """Additional arguments of ``bwa-mem3 mem``"""


class Correct(SnappyModel):
    """``fgumi correct``: the extracted UMIs are corrected against the known UMIs of
    ``depends_on.umis``, for kits with a fixed set of UMIs"""

    max_mismatches: Annotated[int, Field(ge=0)] = 2
    """Mismatches up to which a UMI is corrected to a known one"""

    min_distance: Annotated[int, Field(ge=0)] = 1
    """Mismatches by which the best known UMI must beat the second best"""


class Group(SnappyModel):
    """``fgumi group``: the raw reads of one molecule, by UMI and position"""

    edits: Annotated[int, Field(ge=0)] = 1
    """Mismatches between the UMIs of one molecule; for duplex UMIs over both halves"""

    min_map_q: Annotated[int, Field(ge=0)] = 1
    """Templates with a lower mapping quality are left out"""

    allow_unmapped: bool = False
    """Keep templates whose reads are all unmapped (cfDNA assays)"""

    extra_arguments: list[str] = []
    """Additional arguments of ``fgumi group``"""


class Duplex(SnappyModel):
    """``fgumi duplex`` and the duplex settings of ``fgumi filter``"""

    min_reads: Annotated[list[Annotated[int, Field(ge=0)]], Field(min_length=1, max_length=3)] = [1]
    """Raw reads a consensus read needs: for the duplex read, the AB strand and the BA strand
    (``fgumi duplex --min-reads``); a missing value repeats the last. ``[1]`` keeps every molecule
    for the filter."""

    min_input_base_quality: Annotated[int, Field(ge=0)] = 20
    """Raw bases below this quality are masked before calling"""

    filter_min_reads: Annotated[
        list[Annotated[int, Field(ge=0)]], Field(min_length=1, max_length=3)
    ] = [3, 2, 2]
    """Raw reads a consensus read needs to pass the filter: for the duplex read, the AB strand and
    the BA strand (``fgumi filter --min-reads``)"""

    require_single_strand_agreement: bool = True
    """Mask the bases where the two single-strand consensus reads disagree (``fgumi filter``)"""


class Simplex(SnappyModel):
    """``fgumi simplex`` and the simplex settings of ``fgumi filter``"""

    min_reads: Annotated[int, Field(ge=1)] = 1
    """Raw reads a consensus read needs (``fgumi simplex --min-reads``); 1 keeps every molecule
    for the filter"""

    min_input_base_quality: Annotated[int, Field(ge=0)] = 20
    """Raw bases below this quality are masked before calling"""

    filter_min_reads: Annotated[int, Field(ge=1)] = 3
    """Raw reads a consensus read needs to pass the filter (``fgumi filter --min-reads``)"""


class Filter(SnappyModel):
    """``fgumi filter``, the settings that duplex and simplex share"""

    max_read_error_rate: Annotated[float, Field(ge=0, le=1)] = 0.025
    """Consensus reads whose raw reads disagree with them more often fail"""

    max_base_error_rate: Annotated[float, Field(ge=0, le=1)] = 0.1
    """Consensus bases whose raw bases disagree with them more often are masked"""

    min_base_quality: Annotated[int, Field(ge=0)] = 40
    """Consensus bases below this quality are masked"""

    max_no_call_fraction: Annotated[float, Field(ge=0)] = 0.2
    """Consensus reads with a higher fraction of masked bases fail; a value of 1 or more is a
    count"""

    extra_arguments: list[str] = []
    """Additional arguments of ``fgumi filter``"""


class Realign(SnappyModel):
    """``fgumi realign``: the consensus reads are aligned with bwa-mem3, preferring the window of
    the molecule's raw reads"""

    chunk_size: Annotated[int, Field(gt=0)] = 10_000_000
    """Bases per batch of the realigner's ``bwa-mem3 mem`` (``-K``) and its in-flight budget.
    fgumi's default of 150 Mbp deadlocks on consensus reads with per-base tags."""


class Clip(SnappyModel):
    """``fgumi clip``: hard clipping of read ends after the filter, with the per-base consensus
    tags clipped alike; for example 15 and 8 bases against end-repair artifacts"""

    read_one_five_prime: Annotated[int, Field(ge=0)] = 0
    """Bases clipped from the 5' end of read 1"""

    read_one_three_prime: Annotated[int, Field(ge=0)] = 0
    """Bases clipped from the 3' end of read 1"""

    read_two_five_prime: Annotated[int, Field(ge=0)] = 0
    """Bases clipped from the 5' end of read 2"""

    read_two_three_prime: Annotated[int, Field(ge=0)] = 0
    """Bases clipped from the 3' end of read 2"""

    overlapping_reads: bool = False
    """Clip the overlap of the two reads of a pair, so that a base is counted once"""

    extra_arguments: list[str] = []
    """Additional arguments of ``fgumi clip``"""


class Fgumi(ToolModel):
    """fgumi: UMI extraction, alignment with bwa-mem3, grouping by UMI, consensus calling,
    realignment and filtering"""

    aligner: Aligner = Aligner()

    num_threads: Annotated[int, Field(gt=0)] = 4
    """Threads of each fgumi stage (fastq, zipper, sort, group, the caller, realign, filter,
    clip)"""

    platform: str = "ILLUMINA"
    """Platform (``PL``) of the read groups"""

    compression_level: Annotated[int, Field(ge=0, le=12)] = 6
    """Compression level of the final BAM file"""

    correct: Correct | None = None
    """Correct the UMIs against the known UMIs of ``depends_on.umis``"""

    group: Group = Group()

    duplex: Duplex | None = None
    """The duplex caller and filter, with ``consensus: duplex``"""

    simplex: Simplex | None = None
    """The simplex caller and filter, with ``consensus: simplex``"""

    filter: Filter = Filter()

    realign: Realign = Realign()

    clip: Clip | None = None
    """Clip the read ends of the filtered consensus reads"""

    extra_arguments: list[str] = []
    """Additional arguments of the consensus caller, ``fgumi duplex`` or ``fgumi simplex``"""


class UmiMapping(SnappyStepModel):
    bam_intermediate = True

    depends_on: UmiMappingDependsOn

    tool: Tool

    consensus: Consensus
    """``duplex``: UMIs on both reads, and the reads of both strands of a molecule form one duplex
    consensus read (``fgumi group --strategy paired``, ``fgumi duplex``); ``simplex``: one
    consensus read per UMI family (``--strategy adjacency``, ``fgumi simplex``). The section
    ``fgumi.duplex`` or ``fgumi.simplex`` must be present."""

    fgumi: Fgumi | None = None

    @model_validator(mode="after")
    def validate_consensus_section(self):
        """The section of the selected consensus mode is explicit, like a tool section."""
        if self.fgumi is not None and getattr(self.fgumi, str(self.consensus)) is None:
            raise ValueError(
                f"consensus={self.consensus} requires the section fgumi.{self.consensus}"
            )
        return self
