import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


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

    path_somatic_targeted_seq_cnv_calling: str = ""

    ascat: Ascat | None = None

    @model_validator(mode="after")
    def check_tool_cnv_calling(self):
        if self.tool_cnv_calling == "copywriter" and not self.path_somatic_targeted_seq_cnv_calling:
            raise ValueError(
                "When using 'copywriter' as tool_cnv_calling, "
                "path_somatic_targeted_seq_cnv_calling must be set"
            )
        return self
