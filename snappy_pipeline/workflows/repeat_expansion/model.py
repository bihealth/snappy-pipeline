import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class RepeatExpansionDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class Tool(enum.StrEnum):
    expansionhunter = "expansionhunter"


class RepeatExpansion(SnappyStepModel):
    depends_on: RepeatExpansionDependsOn = Field(default_factory=RepeatExpansionDependsOn)

    tool: Tool = Tool.expansionhunter
    """Tool to use for repeat expansion"""

    repeat_catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    repeat_annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""
