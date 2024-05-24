import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyStepModel, EnumField, SnappyModel


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


class TargetedSeqMeiCalling(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.scramble])]

    scramble: Scramble | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
