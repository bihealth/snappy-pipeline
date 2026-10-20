from typing import Annotated

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    COHORT_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedAnnotatedVariants
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants
from snappy_pipeline.workflows.variant_phasing.model import ExpectedPhasedVariants


class IgvSessionGenerationDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    phased_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "phased"})),
        ExpectedPathSchema(ExpectedPhasedVariants),
    ] = ""
    """Phased variants; set exactly one of the three variant fields"""
    annotated_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "annotated"})),
        COHORT_LEVEL,
        ExpectedPathSchema(ExpectedAnnotatedVariants),
    ] = ""
    """Annotated variants; set exactly one of the three variant fields"""
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        COHORT_LEVEL,
        ExpectedPathSchema(ExpectedGermlineVariants),
    ] = ""
    """Variant calls; set exactly one of the three variant fields"""


class IgvSessionGeneration(SnappyStepModel):
    depends_on: IgvSessionGenerationDependsOn

    @model_validator(mode="after")
    def ensure_at_least_one_dependency_is_specified(self):
        if not any(
            getattr(self.depends_on, path)
            for path in ("phased_variants", "annotated_variants", "variants")
        ):
            raise ValueError("No dependency specified for variant phasing, annotation or calling")
        return self
