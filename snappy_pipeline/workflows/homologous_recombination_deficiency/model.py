import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
)


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


class HomologousRecombinationDeficiencyDependsOn(SnappyModel):
    copy_number: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedSomaticCnvCalls),
    ]


class HomologousRecombinationDeficiency(SnappyStepModel):
    depends_on: HomologousRecombinationDeficiencyDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    scarHRD: ScarHRD | None = None
