import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, RequiredBy, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Mehari, Vep
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
    TOOL,
)
from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    vep = "vep"
    mehari = "mehari"


class ExpectedVariantVcf(SnappyModel):
    """Generic VCF contract consumed by the unified variant_annotation step."""

    vcf: str
    vcf_tbi: str


class ExpectedAnnotatedVariants(SnappyModel):
    """Consumer-driven contract: expected annotated-variant output keys."""

    vcf: str
    vcf_tbi: str


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
        ExpectedPathSchema(ExpectedVariantVcf),
    ]

    reference: Reference

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
        RequiredBy(when="calls_report.igv"),
    ] = ""
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
