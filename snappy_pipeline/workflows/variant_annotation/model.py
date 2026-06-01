import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.annotation import Vep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants


class Tool(enum.StrEnum):
    vep = "vep"


class VepCustom(Vep):
    buffer_size: int = 100000
    num_threads: int = 16

    cache_version: str = "85"
    """The cache version to use.  gnomAD v2 used 85, gnomAD v3.1 uses 101."""

    assembly: str = "GRCh37"
    """The assembly to use.  gnomAD v2 used "GRCh37", gnomAD v3.1 uses "GRCh38"."""

    more_flags: str = "--af_gnomade --af_gnomadg"


class ExpectedAnnotatedGermlineVariants(SnappyModel):
    """Consumer-driven contract: expected output keys from variant_annotation."""

    vcf: str
    vcf_tbi: str


class VariantAnnotationDependsOn(SnappyModel):
    variant_calling: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ] = "variant_calling"


class VariantAnnotation(SnappyStepModel):
    depends_on: VariantAnnotationDependsOn = Field(default_factory=VariantAnnotationDependsOn)

    tool: Annotated[Tool, EnumField(Tool, default=Tool.vep)]

    vep: VepCustom | None = None
