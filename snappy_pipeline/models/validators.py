import pydantic
from pydantic import BaseModel


def validate_ngs_mapping_or_link():
    def validate_depends_on_inputs(instance):
        depends_on = getattr(instance, "depends_on", None)
        if depends_on is None:
            raise ValueError("depends_on configuration is required")

        if not getattr(depends_on, "alignments", "") and not getattr(depends_on, "reads", ""):
            raise ValueError("Either depends_on.alignments or depends_on.reads must be set")
        return instance

    return pydantic.model_validator(mode="after")(validate_depends_on_inputs)


class NgsMappingMixin(BaseModel):
    """
    Validate contract-based upstream mapping/raw-provider dependencies.
    """

    _validate_ngs_mapping_or_link = validate_ngs_mapping_or_link()
