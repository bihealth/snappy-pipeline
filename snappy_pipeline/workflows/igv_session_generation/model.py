from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedAnnotatedVariants
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants
from snappy_pipeline.workflows.variant_phasing.model import ExpectedPhasedVariants


class IgvSessionGenerationDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    variant_phasing: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "phased"})),
        ExpectedPathSchema(ExpectedPhasedVariants),
    ] = ""
    variant_annotation: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "annotated"})),
        ExpectedPathSchema(ExpectedAnnotatedVariants),
    ] = ""
    variant_calling: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ] = ""


class IgvSessionGeneration(SnappyStepModel):
    depends_on: IgvSessionGenerationDependsOn = Field(default_factory=IgvSessionGenerationDependsOn)

    @model_validator(mode="after")
    def ensure_at_least_one_dependency_is_specified(self):
        if not any(
            getattr(self.depends_on, path)
            for path in ("variant_phasing", "variant_annotation", "variant_calling")
        ):
            raise ValueError("No dependency specified for variant phasing, annotation or calling")
        return self
