import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedSomaticCnvCalls(SnappyModel):
    """Consumer-driven contract: expected output keys from somatic CNV caller steps."""

    seqz: str


class Tool(enum.StrEnum):
    scarHRD = "scarHRD"


class GenomeName(enum.StrEnum):
    grch37 = "grch37"
    grch38 = "grch38"
    mouse = "mouse"


class ScarHRD(SnappyModel):
    genome_name: GenomeName = GenomeName.grch37

    chr_prefix: bool = False

    length: int = 50
    """Wiggle track for GC reference file"""


class HomologousRecombinationDeficiencyDependsOn(SnappyModel):
    copy_number: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedSomaticCnvCalls),
    ]

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]


class HomologousRecombinationDeficiency(SnappyStepModel):
    depends_on: HomologousRecombinationDeficiencyDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    scarHRD: ScarHRD | None = None
