import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    samtools = "samtools"


class Samtools(SnappyModel):
    min_X_female: float = 1.5
    """Minimum X ploidy for females"""

    max_Y_female: float = 0.33
    """Maximum Y ploidy for females"""

    min_Y_male: float = 0.5
    """Minimum Y ploidy for males"""

    max_X_male: float = 1.5
    """Maximum X ploidy for males"""


class GuessSex(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.samtools], min_length=1)]

    tool_ngs_mapping: str = "bwa"
    """
    Configuration with read mapper and path to mapping output.
    Will use this for generating coverages over autosomes & sex chromosomes using samtools.
    """

    path_ngs_mapping: str = "../ngs_mapping"

    samtools: Samtools = Samtools()
