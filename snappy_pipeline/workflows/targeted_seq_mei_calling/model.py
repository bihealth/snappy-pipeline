import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    scramble = "scramble"


class Scramble(SnappyModel):
    blast_ref: str
    """path to FASTA reference with BLAST DB (`makeblastdb`)"""

    mei_refs: str | None = None
    """MEI reference file (FASTA), if none provided will use default."""

    n_cluster: int = 5
    """minimum cluster size, depth of soft-clipped reads."""

    mei_score: int = 50
    """minimum MEI alignment score."""

    indel_score: int = 80
    """minimum INDEL alignment score."""

    mei_polya_frac: Annotated[float, Field(ge=0, le=1)] = 0.75
    """minimum fraction of clipped length for calling polyA tail."""


class TargetedSeqMeiCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.scramble], min_length=1)]

    scramble: Scramble | None = None
