import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedVariants,
)


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


class ConcatVariants(SnappyStepModel):
    depends_on: ConcatVariantsDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    bcftools: Bcftools = Field(default_factory=Bcftools)
