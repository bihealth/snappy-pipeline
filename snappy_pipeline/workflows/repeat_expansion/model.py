from snappy_pipeline.models import SnappyStepModel


class RepeatExpansion(SnappyStepModel):
    repeat_catalog: str
    """Repeat expansions definitions - used in ExpansionHunter call"""

    repeat_annotation: str
    """Repeat expansions annotations, e.g., normality range - custom file"""

    path_ngs_mapping: str
    """Path to the ngs_mapping step"""
