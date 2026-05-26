import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.link_in.model import ExpectedLinkedRawFastq
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    optitype = "optitype"
    arcashla = "arcashla"


class Optitype(SnappyModel):
    max_reads: int = 5000
    """5000 is a suggestion by OptiType author"""

    num_mapping_threads: int = 4


class ArcasHla(SnappyModel):
    mapper: str = "star"


class ExpectedHlaTyping(SnappyModel):
    """Consumer-driven contract: expected output keys from hla_typing."""

    done: str


class HlaTypingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    link_in: Annotated[
        str,
        DataSignature(DataType.RAW),
        ExpectedPathSchema(ExpectedLinkedRawFastq),
    ] = ""
    """Optional: name of the ``link_in`` task to use as the preprocessed FASTQ source."""


class HlaTyping(SnappyStepModel):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: HlaTypingDependsOn = Field(default_factory=HlaTypingDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.optitype)]

    optitype: Optitype = Optitype()

    arcashla: ArcasHla = ArcasHla()
