from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class BesenbacherParams(SnappyModel):
    """parameters for Besenbacher quality filter"""

    min_gq: int = 50
    min_dp: int = 10
    max_dp: int = 120
    min_ab: float = 0.20
    max_ab: float = 0.90
    max_ad2: int = 1


class VariantDenovoFiltration(SnappyStepModel):
    path_variant_phasing: str = ""

    path_variant_annotation: str = ""

    path_variant_calling: str = ""

    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    tools_ngs_mapping: list[str] = []
    """defaults to ngs_mapping tool"""

    tools_variant_calling: list[str] = []
    """defaults to variant_annotation tool"""

    info_key_reliable_regions: list[str] = []
    """optional INFO keys with reliable regions"""

    info_key_unreliable_regions: list[str] = []
    """optional INFO keys with unreliable regions"""

    params_besenbacher: BesenbacherParams = BesenbacherParams()

    bad_region_expressions: Annotated[
        list[str], Field(examples=[["'UCSC_CRG_MAPABILITY36 == 1'", "'UCSC_SIMPLE_REPEAT == 1'"]])
    ]

    collect_msdn: bool = True
    """whether or not to collect MSDN (requires GATK HC+UG)"""

    @model_validator(mode="after")
    def ensure_variant_paths_are_configured(self: Self) -> Self:
        assert (
            self.path_variant_phasing or self.path_variant_annotation or self.path_variant_calling
        )
        return self
