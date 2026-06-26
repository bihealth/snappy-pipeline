import enum

from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel


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


class Tools(SnappyModel):
    dna: list[MHCIClassDnaTool | MHCIIClassDnaTool] = []
    rna: list[MHCIClassRnaTool | MHCIIClassRnaTool] = []


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


class HlaLa(SnappyModel):
    mapper: str = "bwa"

    path_graph: str | None = None
    """HLA-LA will use that path when provided. The graphs must be serialized"""

    start: Annotated[int, Field(examples=[28477796, 28510120])] = 0
    """GRCh37: 28477796, GRCh38: 28510120"""

    end: Annotated[int, Field(examples=[33448355, 33480577])] = 0
    """GRCh37: 33448355, GRCh38: 33480577"""

    min_score: float = 0.95


class HlaTyping(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Tools = Tools()

    optitype: Optitype = Optitype()

    arcashla: ArcasHla = ArcasHla()

    hla_la: HlaLa = HlaLa()

    @model_validator(mode="after")
    def ensure_at_least_one_tool(self):
        if len(self.tools.dna) + len(self.tools.rna) == 0:
            raise ValueError("No HLA typing tool requested")
        return self
