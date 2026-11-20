import enum
from typing import Annotated

from pydantic import BaseModel

from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    ANY_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class RenameCombine(enum.StrEnum):
    TUMOR = "tumor"
    GERMLINE = "germline"


class ExpectedSomaticVariant(BaseModel):
    vcf: str
    vcf_tbi: str


class ExpectedGermlineVariant(BaseModel):
    vcf: str
    vcf_tbi: str


class CombineVariantsDependsOn(SnappyModel):
    somatic_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedSomaticVariant),
    ]
    germline_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedGermlineVariant),
        ANY_LEVEL,
    ]
    """The germline calls of the matched normal: its own, or the joint calls of its related group"""

    reference: Reference

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
        RequiredBy(when="calls_report.igv"),
    ] = ""
    """The BAMs of the libraries, for the IGV views of ``calls_report.igv``"""


class CombineVariants(SnappyStepModel, ToolModel):
    depends_on: CombineVariantsDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""
    rename_combined: RenameCombine | None = None
