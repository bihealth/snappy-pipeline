from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class SomaticHlaLohCallingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    hla_types: Annotated[
        str,
        DataSignature(DataType.TABULAR, frozenset({"hla"})),
        ExpectedPathSchema(ExpectedHlaTyping),
    ]


class SomaticHlaLohCalling(SnappyStepModel):
    depends_on: SomaticHlaLohCallingDependsOn

    path_somatic_purity_ploidy: str
