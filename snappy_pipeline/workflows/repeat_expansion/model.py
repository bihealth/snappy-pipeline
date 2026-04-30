from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class RepeatExpansionDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class RepeatExpansion(SnappyStepModel):
    depends_on: RepeatExpansionDependsOn = Field(default_factory=RepeatExpansionDependsOn)

    repeat_catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    repeat_annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""

    """Path to the ngs_mapping step"""
