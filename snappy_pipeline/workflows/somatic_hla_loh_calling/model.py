from typing import Annotated

from snappy_pipeline.models import ConfigPath
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
    """HLA types of the normal library, from an ``hla_typing`` task"""


class SomaticHlaLohCalling(SnappyStepModel):
    depends_on: SomaticHlaLohCallingDependsOn

    path_hla_dat: ConfigPath
    """HLA exon locations, ``hla.dat`` of the OptiType data"""

    path_hla_fasta: ConfigPath
    """HLA reference sequences, ``hla_reference_dna.fasta`` of the OptiType data"""

    path_picard_dir: ConfigPath
    """Directory with the picard-tools 1.x jar files that LOHHLA calls (``--gatkDir``)"""
