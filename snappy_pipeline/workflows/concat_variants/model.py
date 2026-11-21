import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
)
from snappy_pipeline.models.calls_report import CallsReport


class Tool(enum.StrEnum):
    bcftools = "bcftools"


class Bcftools(ToolModel):
    """``bcftools concat --allow-overlaps``; ``extra_arguments`` go to it, such as ``--rm-dups exact``"""


class ConcatVariantsDependsOn(SnappyModel):
    variants: Annotated[
        list[str],
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
    ] = Field(min_length=2)
    """The tasks whose calls are concatenated, for example the small variants and the structural
    variants of each pedigree; they must write the same entities with the same samples"""

    alignments: Annotated[DnaAlignments, RequiredBy(when="calls_report.igv")] = ""
    """The BAMs of the libraries, for the IGV views of ``calls_report.igv``"""

    reference: Annotated[Reference, RequiredBy(when="calls_report.igv")] = ""
    """The genome, for the IGV views of ``calls_report.igv``"""


class ConcatVariants(SnappyStepModel):
    depends_on: ConcatVariantsDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""

    tool: Annotated[Tool, EnumField(Tool)]

    bcftools: Bcftools = Field(default_factory=Bcftools)
