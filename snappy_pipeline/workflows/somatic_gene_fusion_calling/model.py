import enum
import os
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Tool(enum.StrEnum):
    fusioncatcher = "fusioncatcher"
    jaffa = "jaffa"
    arriba = "arriba"
    defuse = "defuse"
    hera = "hera"
    pizzly = "pizzly"
    star_fusion = "star_fusion"


class Fusioncatcher(SnappyModel):
    data_dir: str
    configuration: str = ""
    num_threads: int = 16


class Pizzly(SnappyModel):
    kallisto_index: str
    transcripts_fasta: str
    annotations_gtf: str
    kmer_size: Annotated[int, Field(31, gt=0)]


class Hera(SnappyModel):
    path_index: str
    path_genome: str


class StarFusion(SnappyModel):
    path_ctat_resource_lib: str


class Defuse(SnappyModel):
    path_dataset_directory: str


class Jaffa(SnappyModel):  # TODO
    pass


class Arriba(SnappyModel):
    path_index: str
    """STAR path index (preferably 2.7.10 or later)"""

    blacklist: str = ""
    """provided in the arriba distribution, see /fast/work/groups/cubi/projects/biotools/static_data/app_support/arriba/v2.3.0"""

    known_fusions: str = ""

    tags: str = ""
    """can be set to the same path as known_fusions"""

    structural_variants: str = ""

    protein_domains: str = ""

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

    @model_validator(mode="after")
    def ensure_star_index_files_exist(self):
        full_path = self.path_index
        # a lot of files should be in this dir, justtest these
        for indfile in ("Genome", "SA", "SAindex"):
            expected_path = os.path.join(full_path, indfile)
            if not os.path.exists(expected_path):  # pragma: no cover
                raise ValueError(f"Expected STAR index file {expected_path} does not exist!")
        return self


class SomaticGeneFusionCalling(SnappyStepModel):
    path_link_in: str = ""
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[
        list[Tool],
        EnumField(
            Tool,
            [
                Tool.fusioncatcher,
                Tool.jaffa,
                Tool.arriba,
                Tool.defuse,
                Tool.hera,
                Tool.pizzly,
                Tool.star_fusion,
            ],
        ),
    ]

    fusioncatcher: Fusioncatcher | None = None

    jaffa: Jaffa | None = None

    arriba: Arriba | None = None

    defuse: Defuse | None = None

    hera: Hera | None = None

    pizzly: Pizzly | None = None

    star_fusion: StarFusion | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
