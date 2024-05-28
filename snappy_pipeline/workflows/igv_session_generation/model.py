from typing import Self

from pydantic import model_validator

from snappy_pipeline.models import SnappyStepModel


class IgvSessionGeneration(SnappyStepModel):
    path_ngs_mapping: str

    path_variant_phasing: str = ""

    path_variant_annotation: str = ""

    path_variant_calling: str = ""

    tools_ngs_mapping: list[str] = []
    """defaults to ngs_mapping tool"""

    tools_variant_calling: list[str] = []
    """defaults to variant_annotation tool"""

    @model_validator(mode="after")
    def ensure_at_least_one_path_is_specified(self: Self) -> Self:
        if not any(
            getattr(self, path)
            for path in ("path_variant_phasing", "path_variant_annotation", "path_variant_calling")
        ):
            raise ValueError("No path specified for variant phasing, annotation or calling")
        return self
