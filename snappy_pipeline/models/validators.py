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
    def path_ngs_mapping_or_path_link_in(instance):
        if not instance.path_ngs_mapping and not instance.path_link_in:
            raise ValueError("Either path_ngs_mapping or path_link_in must be set")
        return instance

    return pydantic.model_validator(mode="after")(path_ngs_mapping_or_path_link_in)


class NgsMappingMixin(BaseModel):
    """
    A mixin for validating that not both `path_ngs_mapping` and `path_link_in` are set.
    """

    _validate_ngs_mapping_or_link = validate_ngs_mapping_or_link()
