import pydantic
from pydantic import BaseModel


def validate_ngs_mapping_or_link():
    def validate_depends_on_inputs(instance):
        depends_on = getattr(instance, "depends_on", None)
        if depends_on is None:
            raise ValueError("depends_on configuration is required")

        path_ngs_mapping = getattr(depends_on, "ngs_mapping", "")
        path_link_in = getattr(depends_on, "link_in", "") or getattr(
            depends_on, "adapter_trimming", ""
        )

        if not path_ngs_mapping and not path_link_in:
            raise ValueError(
                "Either depends_on.ngs_mapping or a raw provider dependency "
                "(depends_on.link_in/depends_on.adapter_trimming) must be set"
            )
        return instance

    return pydantic.model_validator(mode="after")(validate_depends_on_inputs)


class NgsMappingMixin(BaseModel):
    """
    Validate contract-based upstream mapping/raw-provider dependencies.
    """

    _validate_ngs_mapping_or_link = validate_ngs_mapping_or_link()
