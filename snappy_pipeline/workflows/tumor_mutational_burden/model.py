from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Kits,
)
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants


class TumorMutationalBurdenDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic"})),
        ExpectedPathSchema(ExpectedSomaticVariants),
    ]

    kits: Kits
    """The TMB counts the variants in the ``targets`` of each tumor library's kit"""


class TumorMutationalBurden(SnappyStepModel):
    depends_on: TumorMutationalBurdenDependsOn

    missense_regex: str = r".*[\\|&]missense_variant[\\|&].*"
    """change if the annotation tool doesn't use 'missense_variant' to indicate missense variant"""

    has_annotation: bool = True
    """Whether the input VCF has been annotated (needed by wrapper for field parsing)."""
