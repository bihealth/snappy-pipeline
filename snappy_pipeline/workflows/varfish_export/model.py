from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Exons,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignmentQc
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants


class ExpectedMehariDatabase(SnappyModel):
    hgnc_tsv: str
    transcript_db: str
    clinvar_db: str
    frequency_db: str


class VarfishExportDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignmentQc),
    ]
    """Used output of ngs_mapping is alignment quality control data"""

    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ]
    """Used output of variant_calling is variant calls"""

    structural_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"})),
    ] = ""
    """Optional ``sv_calling_targeted`` or ``sv_calling_wgs`` task whose calls are exported"""

    reference: Reference

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({"mehari"})),
        ExpectedPathSchema(ExpectedMehariDatabase),
    ]
    """``external_data`` task of type ``database`` tagged ``mehari``: the HGNC table and the Mehari databases"""

    exons: Exons = ""
    """Optional task with a BED file of exons; the export keeps only the small variants near them"""


class VarfishExport(SnappyStepModel):
    """Configuration of the input path enables export from the corresponding pipeline step."""

    depends_on: VarfishExportDependsOn
