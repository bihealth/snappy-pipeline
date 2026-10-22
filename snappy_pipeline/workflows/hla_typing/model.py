import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel, validators
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


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


class Optitype(ToolModel):
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


class ArcasHla(ToolModel):
    population: Population = Population.PRIOR
    min_count: int = 75
    tolerance: float = 1e-7
    max_iterations: int = 1000
    drop_iterations: int | None = None
    drop_threshold: float = 0.1
    zygosity_threshold: float = 0.15
    """Proportion of the major allele's abundance needed to call a gene heterozygous"""
    avg: int | None = None
    """Fragment length mean of single-end reads"""
    std: int | None = None
    """Fragment length standard deviation of single-end reads"""


class ExpectedHlaTyping(SnappyModel):
    """Consumer-driven contract: expected output keys from hla_typing."""

    txt: str
    calls_json: str


class HlaTypingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""
    """Alignments of RNA libraries, read by ``arcashla``"""

    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
    ] = ""
    """FASTQ source: ``data_sets``, an ``external_data`` task or an ``adapter_trimming`` task"""


#: The ``depends_on`` fields each tool reads.
TOOL_DEPENDENCIES = {Tool.optitype: ("reads",), Tool.arcashla: ("alignments",)}


class HlaTyping(SnappyStepModel):
    """Override data set configuration search paths for FASTQ files"""

    depends_on: HlaTypingDependsOn = Field(default_factory=HlaTypingDependsOn)

    tool: Annotated[Tool, EnumField(Tool)]

    optitype: Optitype = Optitype()

    arcashla: ArcasHla = ArcasHla()

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
