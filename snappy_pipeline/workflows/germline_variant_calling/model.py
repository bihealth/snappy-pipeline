import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    gatk4_hc = "gatk4_hc"


class Gatk4HcGvcf(SnappyModel):
    window_length: int = 10000000
    num_threads: int = 16
    allow_seq_dict_incompatibility: bool = False


class GermlineVariantCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.gatk4_hc], min_length=1)]

    ignore_chroms: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]

    gatk4_hc: Gatk4HcGvcf | None = None
