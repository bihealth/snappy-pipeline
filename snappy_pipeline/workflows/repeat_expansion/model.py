import enum

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class RepeatExpansionDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class Tool(enum.StrEnum):
    expansionhunter = "expansionhunter"


class RepeatExpansion(SnappyStepModel):
    depends_on: RepeatExpansionDependsOn = Field(default_factory=RepeatExpansionDependsOn)

    tool: Tool
    """Tool to use for repeat expansion"""

    repeat_catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    repeat_annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""
