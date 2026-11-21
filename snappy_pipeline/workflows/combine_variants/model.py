import enum
from typing import Annotated

from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    ANY_LEVEL,
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
)
from snappy_pipeline.models.calls_report import CallsReport


class RenameCombine(enum.StrEnum):
    TUMOR = "tumor"
    GERMLINE = "germline"


class CombineVariantsDependsOn(SnappyModel):
    somatic_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
    ]
    germline_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
        ANY_LEVEL,
    ]
    """The germline calls of the matched normal: its own, or the joint calls of its related group"""

    reference: Reference

    alignments: Annotated[DnaAlignments, RequiredBy(when="calls_report.igv")] = ""
    """The BAMs of the libraries, for the IGV views of ``calls_report.igv``"""


class CombineVariants(SnappyStepModel, ToolModel):
    depends_on: CombineVariantsDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""
    rename_combined: RenameCombine | None = None
