import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyStepModel


class Tool(enum.StrEnum):
    mantis2_msi = "mantis2_msi"


class SomaticMsiCalling(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mantis2_msi], min_length=1)]

    loci_bed: Annotated[
        str,
        Field(
            examples=[
                "/fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed",
                "/fast/work/groups/cubi/projects/biotools/Mantis/appData/hg38/GRCh38.d1.vd1.all_loci.bed",
            ]
        ),
    ]
