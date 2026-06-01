from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class SomaticHlaLohCallingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    hla_typing: Annotated[
        str,
        DataSignature(DataType.TABULAR, frozenset({"hla"})),
        ExpectedPathSchema(ExpectedHlaTyping),
    ] = "hla_typing"


class SomaticHlaLohCalling(SnappyStepModel):
    depends_on: SomaticHlaLohCallingDependsOn = Field(default_factory=SomaticHlaLohCallingDependsOn)

    path_somatic_purity_ploidy: str
