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


class TargetedSeqMeiCallingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class TargetedSeqMeiCalling(SnappyStepModel):
    depends_on: TargetedSeqMeiCallingDependsOn = Field(
        default_factory=TargetedSeqMeiCallingDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.scramble)]

    scramble: Scramble | None = None
