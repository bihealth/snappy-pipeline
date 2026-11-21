import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, RequiredBy, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
    TOOL,
)
from snappy_pipeline.models.calls_report import CallsReport


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class ExpectedVepDatabase(SnappyModel):
    cache: str
    """The VEP cache directory"""


class ExpectedMehariDatabase(SnappyModel):
    transcript_db: str | list[str]
    """One or more transcript databases"""

    frequency_db: str | None = None

    clinvar_db: str | None = None


class VariantAnnotationDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
    ]

    reference: Reference

    alignments: Annotated[DnaAlignments, RequiredBy(when="calls_report.igv")] = ""
    """The BAMs of the libraries, for the IGV views of ``calls_report.igv``"""

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({TOOL})),
        ExpectedPathSchema.by_tool(vep=ExpectedVepDatabase, mehari=ExpectedMehariDatabase),
    ]
    """``external_data`` task of type ``database`` tagged with the tool: the VEP cache, or the Mehari databases"""


class VariantAnnotation(SnappyStepModel):
    depends_on: VariantAnnotationDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""

    tool: Annotated[Tool, EnumField(Tool)]

    vep: Vep = Field(default_factory=Vep)

    mehari: Mehari = Field(default_factory=Mehari)
