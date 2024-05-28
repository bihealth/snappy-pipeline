import enum

from snappy_pipeline.models import SnappyStepModel, SnappyModel


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


class GeneExpressionQuantification(SnappyStepModel):
    path_link_in: str = ""
    """OPTIONAL Override data set configuration search paths for FASTQ files"""

    tools: list[str] = [
        "strandedness",
        "featurecounts",
        "duplication",
        "dupradar",
        "rnaseqc",
        "stats",
        "salmon",
    ]

    path_ngs_mapping: str

    strand: Strand | int = -1  # TODO: what is this default value of -1?
    """Use 0, 1 or 2 to force unstranded, forward or reverse strand"""

    featurecounts: Featurecounts

    strandedness: Strandedness

    rnaseqc: RnaSeqC

    dupradar: DupRadar

    salmon: Salmon
