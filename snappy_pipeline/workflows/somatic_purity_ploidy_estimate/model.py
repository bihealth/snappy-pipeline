import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    ascat = "ascat"


class Ascat(SnappyModel):
    b_af_loci: str
    """BED file with loci for B allele frequency."""


class SomaticPurityPloidyEstimateDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"


class SomaticPurityPloidyEstimate(SnappyStepModel):
    depends_on: SomaticPurityPloidyEstimateDependsOn = Field(
        default_factory=SomaticPurityPloidyEstimateDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.ascat)]

    ascat: Ascat | None = None
