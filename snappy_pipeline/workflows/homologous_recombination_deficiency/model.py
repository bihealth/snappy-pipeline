import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedSomaticCnvCalls(SnappyModel):
    """Consumer-driven contract: expected output keys from somatic CNV caller steps."""

    vcf: str


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
    cnv_calling: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedSomaticCnvCalls),
    ] = "somatic_targeted_seq_cnv_calling"

    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"somatic"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class HomologousRecombinationDeficiency(SnappyStepModel):
    depends_on: HomologousRecombinationDeficiencyDependsOn = Field(
        default_factory=HomologousRecombinationDeficiencyDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.scarHRD)]

    scarHRD: ScarHRD | None = None
