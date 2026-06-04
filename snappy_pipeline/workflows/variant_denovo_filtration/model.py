from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedAnnotatedVariants
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants
from snappy_pipeline.workflows.variant_phasing.model import ExpectedPhasedVariants


class BesenbacherParams(SnappyModel):
    """parameters for Besenbacher quality filter"""

    min_gq: int = 50
    min_dp: int = 10
    max_dp: int = 120
    min_ab: float = 0.20
    max_ab: float = 0.90
    max_ad2: int = 1


class VariantDenovoFiltrationDependsOn(SnappyModel):
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


class VariantDenovoFiltration(SnappyStepModel):
    depends_on: VariantDenovoFiltrationDependsOn = Field(
        default_factory=VariantDenovoFiltrationDependsOn
    )

    info_key_reliable_regions: list[str] = []
    """optional INFO keys with reliable regions"""

    info_key_unreliable_regions: list[str] = []
    """optional INFO keys with unreliable regions"""

    params_besenbacher: BesenbacherParams = BesenbacherParams()

    bad_region_expressions: Annotated[
        list[str], Field(examples=[["'UCSC_CRG_MAPABILITY36 == 1'", "'UCSC_SIMPLE_REPEAT == 1'"]])
    ] = []

    collect_msdn: bool = True
    """whether or not to collect MSDN (requires GATK HC+UG)"""

    @model_validator(mode="after")
    def ensure_variant_paths_are_configured(self):
        assert (
            self.depends_on.variant_phasing
            or self.depends_on.variant_annotation
            or self.depends_on.variant_calling
        )
        return self
