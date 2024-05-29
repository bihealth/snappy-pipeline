import enum
from typing import Annotated, Self

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, EnumField


class Strand(enum.IntEnum):
    unstranded = 0
    forward = 1
    reverse = 2


class Featurecounts(SnappyModel):
    path_annotation_gtf: str


class Strandedness(SnappyModel):
    path_exon_bed: str
    """needs column 6 with strand info, e.g. CCDS/15/GRCh37/CCDS.bed"""

    threshold: float = 0.85


class RnaSeqC(SnappyModel):
    rnaseqc_path_annotation_gtf: str


class DupRadar(SnappyModel):
    dupradar_path_annotation_gtf: str
    num_threads: int = 8


class Salmon(SnappyModel):
    path_transcript_to_gene: str
    path_index: str
    salmon_params: str = " --gcBias --validateMappings"
    num_threads: int = 16


class Tool(enum.Enum):
    strandedness = "strandedness"
    featurecounts = "featurecounts"
    dupradar = "dupradar"
    rnaseqc = "rnaseqc"
    salmon = "salmon"


class GeneExpressionQuantification(SnappyStepModel):
    path_link_in: str = ""
    """OPTIONAL Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool)] = [Tool.salmon]

    path_ngs_mapping: str

    strand: Strand | int = -1  # TODO: what is this default value of -1?
    """Use 0, 1 or 2 to force unstranded, forward or reverse strand. Use -1 to guess."""

    featurecounts: Featurecounts | None = None

    strandedness: Strandedness | None = None

    rnaseqc: RnaSeqC | None = None

    dupradar: DupRadar | None = None

    salmon: Salmon | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
