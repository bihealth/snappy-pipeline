import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyStepModel
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


class VariantAnnotation(SnappyStepModel):
    path_variant_calling: Annotated[str, Field(examples=["../variant_calling"])]
    """Path to variant calling"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep], min_length=1)]

    vep: VepCustom | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
