import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.annotation import Vep


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


class VariantAnnotationDependsOn(SnappyModel):
    variant_calling: str = "variant_calling"


class VariantAnnotation(SnappyStepModel, validators.ToolsMixin):
    depends_on: VariantAnnotationDependsOn = Field(default_factory=VariantAnnotationDependsOn)

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    vep: VepCustom | None = None
