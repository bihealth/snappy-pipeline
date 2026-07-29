import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.link_in.model import ExpectedLinkedRawFastq
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class MHCIClassDnaTool(enum.StrEnum):
    optitype = "optitype"
    hla_la = "hla_la"


class MHCIClassRnaTool(enum.StrEnum):
    optitype = "optitype"
    arcashla = "arcashla"


class MHCIIClassDnaTool(enum.StrEnum):
    hla_la = "hla_la"


class MHCIIClassRnaTool(enum.StrEnum):
    arcashla = "arcashla"


class Tool(enum.StrEnum):
    optitype = "optitype"
    arcashla = "arcashla"


class YaraSensitivity(enum.StrEnum):
    FULL = "full"
    HIGH = "high"
    LOW = "low"


class Yara(SnappyModel):
    error_rate: int = 5
    strata_rate: int = 0
    sensitivity: YaraSensitivity = YaraSensitivity.HIGH


class Optitype(SnappyModel):
    yara_mapper: Yara = Yara()
    max_reads: int = 5000
    """5000 is a suggestion by OptiType author"""
    num_mapping_threads: int = 4
    use_discordant: bool = False


class Population(enum.StrEnum):
    PRIOR = "prior"
    ASIAN_PACIFIC_ISLANDER = "asian_pacific_islander"
    BLACK = "black"
    CAUCASIAN = "caucasian"
    HISPANIC = "hispanic"
    NATIVE_AMERICAN = "native_american"


class ArcasHla(SnappyModel):
    mapper: str = "star"
    population: Population = Population.PRIOR
    min_count: int = 75
    tolerance: float = 1e-7
    max_iterations: int = 1000
    drop_iterations: int | None = None
    drop_threshold: float = 0.1
    zygocity_threshold: float = 0.15
    avg: int | None = None
    std: int | None = None


class ExpectedHlaTyping(SnappyModel):
    """Consumer-driven contract: expected output keys from hla_typing."""

    txt: str
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
