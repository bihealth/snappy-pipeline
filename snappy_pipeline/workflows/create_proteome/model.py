from typing import Annotated

from pydantic import BaseModel, Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema


class ExpectedVariantVcf(BaseModel):
    vcf: str
    vcf_tbi: str


class CreateProteomeDependsOn(SnappyModel):
    variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ]


class CreateProteome(SnappyStepModel):
    depends_on: CreateProteomeDependsOn

    add_reference: bool = False
    """Should the process also compute proteome based on unmutated sequence"""

    path_proteome: Annotated[
        str | None,
        Field(examples=["gencode.v43.pc_translations.fa.gz", "Homo_sapiens.GRCh38.pep.all.fa.gz"]),
    ] = None
    """Path to curated unmutated proteome"""
