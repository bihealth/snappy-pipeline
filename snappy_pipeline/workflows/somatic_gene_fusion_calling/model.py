import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.workflows.abstract.protocol import (
    TOOL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Features,
    Reference,
)
from snappy_pipeline.workflows.reference_index.model import ExpectedIndex


class Tool(enum.StrEnum):
    jaffa = "jaffa"
    arriba = "arriba"
    defuse = "defuse"
    hera = "hera"
    pizzly = "pizzly"
    star_fusion = "star_fusion"


class Pizzly(SnappyModel):
    kallisto_index: ConfigPath
    transcripts_fasta: ConfigPath
    annotations_gtf: ConfigPath
    kmer_size: Annotated[int, Field(31, gt=0)]


class Hera(SnappyModel):
    path_index: ConfigPath
    path_genome: ConfigPath


class StarFusion(SnappyModel):
    pass


class Defuse(SnappyModel):
    path_dataset_directory: ConfigPath


class Jaffa(SnappyModel):
    path_reference_files: ConfigPath
    """Directory with the JAFFA reference files (``JAFFA_REF_BASE``)"""


class Arriba(SnappyModel):
    num_threads: int = 8

    trim_adapters: bool = False

    num_threads_trimming: int = 2

    star_parameters: list[str] = [
        " --outFilterMultimapNmax 50",
        " --peOverlapNbasesMin 10",
        " --alignSplicedMateMapLminOverLmate 0.5",
        " --alignSJstitchMismatchNmax 5 -1 5 5",
        " --chimSegmentMin 10",
        " --chimOutType WithinBAM HardClip",
        " --chimJunctionOverhangMin 10",
        " --chimScoreDropMax 30",
        " --chimScoreJunctionNonGTAG 0",
        " --chimScoreSeparation 1",
        " --chimSegmentReadGapMax 3",
        " --chimMultimapNmax 50",
    ]


class ExpectedArribaDatabase(SnappyModel):
    """Consumer-driven contract: the optional files of the arriba release."""

    blacklist: str = ""
    """Provided in the arriba distribution (``database/`` of the release)"""

    known_fusions: str = ""

    tags: str = ""
    """Can be set to the same path as ``known_fusions``"""

    structural_variants: str = ""

    protein_domains: str = ""


class ExpectedStarFusionDatabase(SnappyModel):
    """Consumer-driven contract: the STAR-Fusion database."""

    ctat_resource_lib: str
    """The CTAT genome resource library"""


class SomaticGeneFusionCallingDependsOn(SnappyModel):
    reads: Annotated[
        str,
        DataSignature(DataType.RAW),
    ]
    """
    FASTQ source: an ``external_data`` task, an ``adapter_trimming`` task,
    or ``data_sets`` to search the data sets' search paths.
    """

    reference: Reference = ""

    features: Features = ""

    index: Annotated[
        str,
        DataSignature(DataType.INDEX, frozenset({"star", "rna"})),
        ExpectedPathSchema(ExpectedIndex),
    ] = ""
    """Needed for ``tool: arriba``: the STAR index (preferably 2.7.10 or later)"""

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({TOOL})),
        ExpectedPathSchema.by_tool(
            arriba=ExpectedArribaDatabase, star_fusion=ExpectedStarFusionDatabase
        ),
    ] = ""
    """The database of the task's tool: optional for ``tool: arriba``, needed for
    ``tool: star_fusion``"""


TOOL_DEPENDENCIES = {
    Tool.arriba: ("reference", "features", "index"),
    Tool.star_fusion: ("database",),
}


class SomaticGeneFusionCalling(SnappyStepModel):
    depends_on: SomaticGeneFusionCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    jaffa: Jaffa | None = None

    arriba: Arriba | None = None

    defuse: Defuse | None = None

    hera: Hera | None = None

    pizzly: Pizzly | None = None

    star_fusion: StarFusion | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
