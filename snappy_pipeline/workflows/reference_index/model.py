from enum import StrEnum
from typing import Annotated

from pydantic import (
    BaseModel,
    model_validator,
)

from snappy_pipeline.models import (
    SnappyModel,
    SnappyStepModel,
    ToolModel,
)
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Features,
)
from snappy_pipeline.workflows.reference_download.model import (
    ExpectedReferenceDownloadFiles,
    Molecule,
)


class ExpectedIndex(BaseModel):
    """Consumer-driven contract: the index of one mapper.

    The file prefix for bwa and bwa-mem2, the ``.mmi`` file for minimap2, the directory for STAR.
    """

    index: str


class Tool(StrEnum):
    bwa = "bwa"
    bwa_mem2 = "bwa_mem2"
    minimap2 = "minimap2"
    star = "star"


class Bwa(ToolModel):
    algorithm: str = "bwtsw"


class BwaMem2(ToolModel):
    pass


class Minimap2(ToolModel):
    pass


class Star(ToolModel):
    pass


class ReferenceIndexDependsOn(SnappyModel):
    reference: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"reference", ("dna", "rna")})),
        ExpectedPathSchema(ExpectedReferenceDownloadFiles),
    ]
    """Reference FASTA: a ``reference_download`` task, or an ``external_data`` task with the key ``fasta``."""

    features: Features = ""
    """Gene annotation for the STAR index."""


class ReferenceIndex(SnappyStepModel):
    depends_on: ReferenceIndexDependsOn

    tool: Tool
    """Index family to build for this task (one tool per task)."""

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
