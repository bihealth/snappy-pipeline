import enum
from typing import Annotated

from pydantic import model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    ascat = "ascat"


class Ascat(SnappyModel):
    b_af_loci: str
    """BED file with loci for B allele frequency."""


class SomaticPurityPloidyEstimate(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.ascat], min_length=1)]

    tool_cnv_calling: str = "cnvetti"
    """When set to 'copywriter', will trigger 'somatic_targeted_seq_cnv_calling'"""

    tool_ngs_mapping: str = "bwa"
    """
        Configuration with read mapper and path to mapping output.
        Will use this for generating a pileup using samtools
        for obtaining the b allele fraction and computing coverage.
    """

    path_ngs_mapping: str = "../ngs_mapping"

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
