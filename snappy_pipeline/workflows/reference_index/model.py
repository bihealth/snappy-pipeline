from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, model_validator

from pydantic import Field

from snappy_pipeline.models import ResolvablePath, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.reference_download.model import (
    ExpectedReferenceDownloadFiles,
    Molecule,
)


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


class ReferenceIndexDependsOn(SnappyModel):
    reference_download: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"reference", ("dna", "rna")})),
        ExpectedPathSchema(ExpectedReferenceDownloadFiles),
    ] = ""


class ReferenceIndex(SnappyStepModel):
    depends_on: ReferenceIndexDependsOn = Field(default_factory=ReferenceIndexDependsOn)

    tool: Tool = Tool.bwa
    """Index family to build for this task (one tool per task)."""

    path_reference: ResolvablePath = ""
    """Optional FASTA path override. Falls back to static_data_config.reference.path when empty."""

    reference_molecule: Molecule = Molecule.dna
    """Molecule class of the input reference used for index generation."""

    bwa: Bwa = Bwa()
    bwa_mem2: BwaMem2 = BwaMem2()
    minimap2: Minimap2 = Minimap2()
    star: Star = Star()

    @model_validator(mode="after")
    def validate_tool_vs_reference_molecule(self):
        if self.tool == Tool.star and self.reference_molecule != Molecule.rna:
            raise ValueError("tool='star' requires reference_molecule='rna'")
        if self.tool != Tool.star and self.reference_molecule != Molecule.dna:
            raise ValueError("tool in {bwa,bwa_mem2,minimap2} requires reference_molecule='dna'")
        return self
