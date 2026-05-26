import pydantic
from pydantic import BaseModel


def validate_tool():
    def ensure_tool_is_configured(instance):
        for tool in [instance.tool]:
            if not getattr(instance, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return instance

    return pydantic.model_validator(mode="after")(ensure_tool_is_configured)


class ToolMixin(BaseModel):
    """
    A mixin for validating that the tool in `self.tool`
    has an accompanying configuration field in the model.
    """

    _validate_tool = validate_tool()


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
