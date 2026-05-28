from enum import StrEnum

from pydantic import BaseModel

from snappy_pipeline.models import ResolvablePath, SnappyModel, SnappyStepModel


class ExpectedReferenceIndexFiles(BaseModel):
    """Consumer-driven contract for reference index artifact paths."""

    bwa_index_prefix: str
    bwa_mem2_index_prefix: str
    minimap2_index: str
    star_index_dir: str
    reference_fai: str
    reference_dict: str
    reference_genome: str


class Tool(StrEnum):
    bwa = "bwa"
    bwa_mem2 = "bwa_mem2"
    minimap2 = "minimap2"
    star = "star"


class Bwa(SnappyModel):
    algorithm: str = "bwtsw"


class BwaMem2(SnappyModel):
    extra_args: list[str] = []


class Minimap2(SnappyModel):
    extra_args: list[str] = []


class Star(SnappyModel):
    extra_args: list[str] = []


class ReferenceIndex(SnappyStepModel):
    tool: Tool = Tool.bwa
    """Index family to build for this task (one tool per task)."""

    path_reference: ResolvablePath = ""
    """Optional FASTA path override. Falls back to static_data_config.reference.path when empty."""

    bwa: Bwa = Bwa()
    bwa_mem2: BwaMem2 = BwaMem2()
    minimap2: Minimap2 = Minimap2()
    star: Star = Star()
