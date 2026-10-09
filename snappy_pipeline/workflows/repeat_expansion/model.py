import enum
from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class RepeatExpansionDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference


class Tool(enum.StrEnum):
    expansionhunter = "expansionhunter"


class RepeatExpansion(SnappyStepModel):
    depends_on: RepeatExpansionDependsOn

    tool: Tool
    """Tool to use for repeat expansion"""

    repeat_catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    repeat_annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""
