import enum
from typing import Annotated, Self

from pydantic import model_validator

from snappy_pipeline.models import SnappyStepModel, EnumField, SnappyModel


class Tool(enum.StrEnum):
    ascat = "ascat"


class Ascat(SnappyModel):
    b_af_loci: str
    """BED file with loci for B allele frequency."""


class SomaticPurityPloidyEstimate(SnappyStepModel):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.ascat])]

    tool_cnv_calling: str = "cnvetti"
    """When set to 'copywriter', will trigger 'somatic_targeted_seq_cnv_calling'"""

    tool_ngs_mapping: str = "bwa"
    """
        Configuration with read mapper and path to mapping output.
        Will use this for generating a pileup using samtools
        for obtaining the b allele fraction and computing coverage.
    """

    path_ngs_mapping: str

    ascat: Ascat | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
