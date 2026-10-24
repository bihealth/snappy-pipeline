import enum
from typing import Annotated

from pydantic import model_validator

from snappy_pipeline.models import ConfigPath, SnappyModel, SnappyStepModel, ToolModel, validators
from snappy_pipeline.workflows.abstract.protocol import (
    CommonVariants,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Kits,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    gatk4_bqsr = "gatk4_bqsr"
    agent = "agent"


class AlignmentProcessingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    """The BAM files to process: an ``ngs_mapping`` task or another ``alignment_processing`` task"""

    reference: Reference

    common_variants: CommonVariants = ""
    """Known variant sites, which BQSR does not count as sequencing errors (``gatk4_bqsr``)"""

    kits: Kits = ""
    """The enrichment kits; AGeNT reads the ``baits`` of each library's kit (``agent``)"""


class Gatk4Bqsr(ToolModel):
    """Base quality score recalibration with GATK ``BaseRecalibrator`` and ``ApplyBQSR``;
    ``extra_arguments`` go to ``BaseRecalibrator``"""


class AgentConsensusMode(enum.StrEnum):
    SINGLE = "SINGLE"
    HYBRID = "HYBRID"
    DUPLEX = "DUPLEX"


class Agent(ToolModel):
    """Duplicate marking with the molecular barcodes (MBCs) of Agilent's AGeNT (LocatIt), from
    reads that AGeNT trimmed (``adapter_trimming`` tool ``agent``) and that the mapper wrote with
    their barcodes (``extra_arguments: ["-C"]`` of the BWA mapper)"""

    path: ConfigPath
    """The jar file of AGeNT's duplicate marking; AGeNT is not open source and not in bioconda"""

    consensus_mode: AgentConsensusMode
    """``-c``"""

    input_filter_args: list[str] = []
    """Arguments after ``-f``, for example ``-mm 13`` (minimum base quality), ``-mr 13`` (minimum
    barcode base quality) and ``-mq 30`` (minimum mapping quality)"""

    consensus_filter_args: list[str] = []
    """Arguments after ``-F``"""

    java_options: str = "-Xmx72000M"
    """Options of the Java virtual machine; the heap is 75% of the job's memory"""

    extra_arguments: list[str] = []
    """For example ``-d 1`` (maximum number of barcode mismatches)"""


TOOL_DEPENDENCIES = {Tool.gatk4_bqsr: ("common_variants",), Tool.agent: ("kits",)}


class AlignmentProcessing(SnappyStepModel):
    bam_intermediate = True

    depends_on: AlignmentProcessingDependsOn

    tool: Tool

    gatk4_bqsr: Gatk4Bqsr | None = None
    agent: Agent | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
