import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from models import SnappyStepModel, EnumField
from models.annotation import Vep


class Tool(enum.Enum):
    vep = "vep"


class VepCustom(Vep):
    buffer_size = 100000
    num_threads = 16

    cache_version = "85"
    """The cache version to use.  gnomAD v2 used 85, gnomAD v3.1 uses 101."""

    assembly = "GRCh37"
    """The assembly to use.  gnomAD v2 used "GRCh37", gnomAD v3.1 uses "GRCh38"."""

    more_flags: str = "--af_gnomade --af_gnomadg"


class VariantAnnotation(SnappyStepModel):
    path_variant_calling: Annotated[str, Field(examples=["../variant_calling"])]
    """Path to variant calling"""

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.vep])]

    vep: VepCustom | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
