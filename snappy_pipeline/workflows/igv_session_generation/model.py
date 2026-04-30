from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class IgvSessionGenerationDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    variant_phasing: str = ""
    variant_annotation: str = ""
    variant_calling: str = ""


class IgvSessionGeneration(SnappyStepModel):
    depends_on: IgvSessionGenerationDependsOn = Field(default_factory=IgvSessionGenerationDependsOn)

    tools_ngs_mapping: list[str] = []
    """defaults to ngs_mapping tool"""

    tools_variant_calling: list[str] = []
    """defaults to variant_annotation tool"""

    @model_validator(mode="after")
    def ensure_at_least_one_dependency_is_specified(self):
        if not any(
            getattr(self.depends_on, path)
            for path in ("variant_phasing", "variant_annotation", "variant_calling")
        ):
            raise ValueError("No dependency specified for variant phasing, annotation or calling")
        return self
