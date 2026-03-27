import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


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
    path_index: str
    salmon_params: str = " --gcBias --validateMappings"
    num_threads: int = 16


class Duplication(SnappyModel):
    pass


class Stats(SnappyModel):
    pass


class Tool(enum.Enum):
    strandedness = "strandedness"
    featurecounts = "featurecounts"
    dupradar = "dupradar"
    duplication = "duplication"
    rnaseqc = "rnaseqc"
    salmon = "salmon"
    stats = "stats"


class GeneExpressionQuantification(
    SnappyStepModel, validators.NgsMappingMixin, validators.ToolsMixin
):
    path_ngs_mapping: str = "../ngs_mapping"

    path_link_in: str = ""
    """OPTIONAL Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, min_length=1)] = [Tool.salmon]

    strand: Strand | int = -1  # TODO: what is this default value of -1?
    """Use 0, 1 or 2 to force unstranded, forward or reverse strand. Use -1 to guess."""

    featurecounts: Featurecounts | None = None

    strandedness: Strandedness | None = None

    rnaseqc: RnaSeqC | None = None

    dupradar: DupRadar | None = None

    duplication: Duplication | None = None

    stats: Stats | None = None

    salmon: Salmon | None = None
