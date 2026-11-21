import enum
from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    TOOL,
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    Reference,
)


class ExpectedExpansionHunterDatabase(SnappyModel):
    """Consumer-driven contract: the files of ExpansionHunter."""

    catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""


class RepeatExpansionDependsOn(SnappyModel):
    alignments: DnaAlignments

    reference: Reference

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({TOOL})),
        ExpectedPathSchema(ExpectedExpansionHunterDatabase),
    ]
    """The repeat catalog and the annotation of the repeats"""


class Tool(enum.StrEnum):
    expansionhunter = "expansionhunter"


class RepeatExpansion(SnappyStepModel, ToolModel):
    depends_on: RepeatExpansionDependsOn

    tool: Tool
    """Tool to use for repeat expansion"""
