import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedRegions,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    mantis_msi2 = "mantis_msi2"


class SomaticMsiCallingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    microsatellites: Annotated[
        str,
        DataSignature(DataType.REGIONS, frozenset({"microsatellites"})),
        ExpectedPathSchema(ExpectedRegions),
    ]
    """Task with the BED file of the microsatellite loci, such as ``appData/hg38/GRCh38.d1.vd1.all_loci.bed`` of MANTIS"""


class SomaticMsiCalling(SnappyStepModel, ToolModel):
    depends_on: SomaticMsiCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    min_read_length: int = 35
    """Minimum read length (``mantis-msi2 --min-read-length``)"""

    min_read_quality: float = 20.0
    """Minimum mean base quality of a read (``--min-read-quality``)"""

    min_locus_quality: float = 25.0
    """Minimum mean base quality of the reads at a locus (``--min-locus-quality``)"""

    min_locus_coverage: int = 20
    """Minimum coverage of a locus in tumor and normal (``--min-locus-coverage``)"""

    min_repeat_reads: int = 1
    """Minimum number of reads supporting a repeat length (``--min-repeat-reads``)"""
