from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants


class TumorMutationalBurdenDependsOn(SnappyModel):
    somatic_variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic"})),
        ExpectedPathSchema(ExpectedSomaticVariants),
    ] = "somatic_variant"


class TumorMutationalBurden(SnappyStepModel):
    depends_on: TumorMutationalBurdenDependsOn = Field(
        default_factory=TumorMutationalBurdenDependsOn
    )

    target_regions: str
    """Path to target_regions file (bed format)"""

    missense_regex: str = r".*[\\|&]missense_variant[\\|&].*"
    """change if the annotation tool doesn't use 'missense_variant' to indicate missense variant"""

    has_annotation: bool = True
    """Whether the input VCF has been annotated (needed by wrapper for field parsing)."""
