import enum
import json
import os
import re
import types
import typing
from enum import Enum
from inspect import isclass
from io import StringIO
from typing import Annotated

import ruamel
import typing_extensions
from annotated_types import Predicate
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, model_validator
from pydantic_core import PydanticUndefined
from ruamel.yaml import YAML


class _PathMarker:
    """Marker for path fields that should be resolved. Used with Annotated."""

    def __init__(self, check_exists: bool = True):
        self.check_exists = check_exists

    def __repr__(self):
        return f"_PathMarker(check_exists={self.check_exists})"


# Public path annotation types for use in models
ResolvablePath = Annotated[str, _PathMarker(check_exists=True)]
"""Annotated type for file paths that must exist after resolution."""

ResolvablePathPrefix = Annotated[str, _PathMarker(check_exists=False)]
"""Annotated type for index prefixes or paths where sidecar files are checked separately."""

ResolvablePathList = Annotated[list[str], _PathMarker(check_exists=True)]
"""Annotated type for lists of file paths that must exist after resolution."""


def resolve_relative_path(path_value: str, info: ValidationInfo, check_exists: bool = True) -> str:
    """
    Resolve relative paths using config_lookup_paths from validation context.

    Args:
        path_value: The path to resolve (string)
        info: ValidationInfo from Pydantic
        check_exists: If True, raises error if relative path cannot be found in any lookup base.
                     If False, returns the path resolved against the first lookup base
                     (useful for index prefixes that don't exist as files themselves).

    Returns:
        Absolute path. The original value is returned if:
          - it's already absolute
          - it's empty
          - context doesn't provide config_lookup_paths

    Raises:
        ValueError if check_exists=True and relative path cannot be found in any lookup base.
    """
    if not path_value or os.path.isabs(path_value):
        return path_value

    config_lookup_paths = info.context.get("config_lookup_paths") if info.context else None
    if not config_lookup_paths:
        return path_value

    # Try to find the path in any lookup base
    for base_path in config_lookup_paths:
        candidate = os.path.abspath(os.path.join(base_path, path_value))
        if os.path.exists(candidate):
            return candidate

    # If check_exists is False, return the path resolved against the first lookup base
    if not check_exists:
        return os.path.abspath(os.path.join(config_lookup_paths[0], path_value))

    # Relative path not found and check_exists=True
    raise ValueError(
        f"relative path '{path_value}' not found in lookup bases: {config_lookup_paths}"
    )


def enum_options(enum: Enum) -> list[tuple[str, typing.Any]]:
    """Returns a list of tuples containing the name and value of each enum member."""
    return [(e.name, e.value) for e in enum]


def EnumField(enum: type[Enum], default: typing.Any = PydanticUndefined, *args, **kwargs):
    """
    An extension of pydantic's `Field` that adds 'options' to the json_schema_extra field,
    containing the available options of the specified enum.
    """
    extra = kwargs.get("json_schema_extra", {})
    extra.update(dict(options=enum_options(enum)))
    kwargs["json_schema_extra"] = extra
    return Field(default, *args, **kwargs)


size_string_regexp = re.compile(r"[. 0-9]+([KMGTP])")
SizeString = Annotated[str, Predicate(lambda s: size_string_regexp.match(s) is not None)]
"""A string representing a size, e.g. '1G' for 1 gigabyte."""


class KeepTmpdir(enum.StrEnum):
    """Whether to keep the temporary directory after the job has finished."""

    always = "always"
    never = "never"
    onerror = "onerror"


class SnappyModel(BaseModel):
    """
    Base class for all snappy models.
    By default, extra fields are forbidden, attribute docstrings are used for field descriptions,
    enum member values instead of names are used, and default values are validated (because
    validation can potentially modify the values of fields with default values)

    Path fields can be marked with ResolvablePath or ResolvablePathPrefix annotations
    to automatically resolve relative paths during validation
    """

    model_config = ConfigDict(
        extra="forbid",
        use_attribute_docstrings=True,
        use_enum_values=True,
        validate_default=True,
        populate_by_name=True,
    )

    def get(self, key: str, default: typing.Any = None) -> typing.Any:
        """
        Return the value of the field with the given key, or the default value if it doesn't exist.
        Simply delegates to getattr.
        """
        return getattr(self, key, default)

    def __getitem__(self, item: str) -> typing.Any:
        """
        Return the value of the field with the given key.
        Raise an AttributeError if the field doesn't exist.
        """
        return getattr(self, item)

    def keys(self):
        """Return a list of field names."""
        return self.model_fields.keys()


class RelationshipDefinition(SnappyModel):
    """Definition of a relationship between rows in the library DataFrame.

    A relationship adds a new column to the DataFrame whose value is looked
    up from a *related* row.  The ``via`` column is the join key (e.g.
    ``"donor_name"``), and ``target`` is a pandas ``DataFrame.query()``
    expression applied to the related rows to select the correct match.

    **Examples**

    .. code-block:: yaml

        # Find the matched normal DNA library for a tumor sample
        matched_normal_lib:
          via: donor_name
          target: "role == 'normal' and extraction_type == 'dna'"
          column: matched_normal_lib

        # Find the index/proband library for any library in the same cohort
        index_lib:
          via: cohort_name
          target: "role == 'index'"
          column: index_lib

        # Find all tumor libraries in the same donor (one-to-many)
        donor_tumor_libs:
          via: donor_name
          target: "role == 'tumor'"
          column: donor_tumor_libs
          many: true
    """

    via: str
    """Column name in the library DataFrame to use as the join key.
    The related row must have the *same* value in this column as the
    source row."""

    target: str
    """Pandas ``DataFrame.query()`` expression applied to the related
    rows (those sharing the same ``via`` value) to select the desired
    match."""

    column: str | None = None
    """Name of the new column to add.  Defaults to the relationship key
    name in the ``relationships`` dict."""

    many: bool = False
    """If ``True``, the relationship may match multiple related rows.
    The column value will be a ``list[str]`` of all matching library
    names.  If ``False`` (default), only the first match is kept and
    the column value is a single ``str`` (or ``""`` if no match)."""


# This exists to distinguish workflow step_config models from other snappy specific models
# It also provides a default_config_yaml_string method that includes the step_config section
# by default.
class SnappyStepModel(SnappyModel, object):
    """A base class for all workflow step configuration models.

    Every step gets ``library_selection``, ``group_by``, and ``relationships``
    fields.  Steps that do not use these fields leave them at their defaults
    (``None`` / empty).
    """

    library_selection: str | None = None
    """Optional pandas ``DataFrame.query()`` expression to choose which
    libraries this task processes.

    The expression is evaluated against the tidy library DataFrame produced
    by :func:`~snappy_pipeline.workflows.abstract.build_library_dataframe`.
    See the ``_LIBRARY_SELECTION_DEFAULTS`` dict and column reference in
    ``build_library_dataframe`` for supported columns."""

    group_by: str | None = None
    """Controls output granularity.

    * ``"cohort"`` -- one output per cohort / pedigree / donor.
    * ``None`` (default) -- one output per library."""

    relationships: dict[str, RelationshipDefinition] | None = None
    """Named relationships that add derived columns to the library DataFrame
    *before* ``library_selection`` is applied.

    Each key is the column name; the value is a
    :class:`RelationshipDefinition`."""

    @model_validator(mode="after")
    def validate_selected_tool_config(self):
        model_fields = type(self).model_fields
        if "tool" not in model_fields:
            return self

        selected = getattr(self, "tool", None)
        if selected is None:
            return self
        if isinstance(selected, enum.Enum):
            selected = selected.value
        if not isinstance(selected, str) or not selected:
            return self

        if selected in model_fields and getattr(self, selected, None) is None:
            raise ValueError(f"tool={selected} requires explicit '{selected}' config section")

        return self

    @classmethod
    def default_config_yaml_string(
        cls, comment_optional: bool = True, with_step_config: bool = True
    ):
        config_str = default_config_yaml_string(cls, comment_optional)
        if with_step_config:
            config_str = (
                (" " * INDENTATION * 2 + line) for line in config_str.splitlines(keepends=True)
            )

            def camel_to_snake(name: str) -> str:
                name = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
                return re.sub("([a-z0-9])([A-Z])", r"\1_\2", name).lower()

            name = camel_to_snake(cls.__name__)
            config_str = f"step_config:\n{' ' * INDENTATION}{name}:\n{''.join(config_str)}"

        return config_str

    def model_dump_yaml(self, **kwargs) -> str:
        cfg, _ = _model_to_commented_yaml(self, **kwargs)
        return _dump_yaml(cfg)


# Classes and functions for generating an annotated default configuration YAML
# from a snappy pydantic model
# list all optional fields commented out, list enum options, mark required fields, show descriptions

INDENTATION = 2  # Indentation used for YAML


def default_config_yaml_string(model: type[SnappyStepModel], comment_optional: bool = False) -> str:
    return _dump_commented_yaml(model, comment_optional)


def _check_model_class(annotation, clazz: type = BaseModel) -> bool:
    """Checks whether the given annotation is a class and is a subclass of `clazz`"""
    try:
        is_model_class = isclass(annotation) and issubclass(annotation, clazz)
    except TypeError:
        # for some reason, types.GenericAlias doesn't count as a class?
        is_model_class = hasattr(annotation, "model_fields")
    return is_model_class


def _annotate_model(
    config_model: type[BaseModel],
    comment_map: ruamel.yaml.CommentedMap,
    max_column: int = 80,
    level: int = 0,
    indent: int = INDENTATION,
    path: list[str] = [],
):
    """
    Annotates a given pydantic model with comments
    based on whether a field is required or not, as well as listing options for enum valued fields
    """
    for key, field in config_model.model_fields.items():
        annotation = field.annotation
        is_union = _is_union_type(annotation)
        allows_none = types.NoneType in typing_extensions.get_args(annotation)
        is_optional = (is_union and allows_none) or not field.is_required()

        if is_optional:
            tags = []
        else:
            tags = ["REQUIRED"]

        comment = tags

        if field.examples:
            comment.append(f"Examples: {', '.join(map(str, field.examples))}")

        options = (getattr(field, "json_schema_extra") or {}).get("options", [])
        if _check_model_class(annotation, enum.Enum):
            options = enum_options(annotation)

        if options:

            def option_str(option: tuple[str, typing.Any]) -> str:
                name, value = option
                if name.upper() != str(value).upper().replace("-", "_"):
                    return f"{repr(value)} ({name})"
                else:
                    return f"{repr(value)}"

            comment.append(f"Options: {', '.join(map(option_str, options))}")

        comment = "; ".join(comment)

        if comment_map is not None:
            if comment:
                comment_map.yaml_add_eol_comment(comment, key, column=max_column)

            if field.description:
                (
                    comment_map.yaml_set_comment_before_after_key(
                        key,
                        indent=indent * level,
                        before="\n" + field.description,
                        after=None,
                    )
                )

        if _check_model_class(annotation):
            _annotate_model(
                annotation,
                comment_map[key],
                level=level + 1,
                max_column=max_column,
                path=path + [key],
            )
        elif _is_union_type(annotation):
            if len(args := typing.get_args(annotation)) == 2 and types.NoneType in args:
                sub_model = next(filter(lambda s: s is not types.NoneType, args))
                if _check_model_class(sub_model):
                    # when the default is set to None but the annotation inherits from basemodel,
                    # make sure to generate those entries as well
                    if not comment_map[key]:
                        sub_model_yaml, m = _model_to_commented_yaml(
                            _placeholder_model_instance(sub_model)
                        )
                        max_column = max(m + 2, max_column)
                        comment_map[key] = sub_model_yaml
                    _annotate_model(
                        sub_model,
                        comment_map[key],
                        level=level + 1,
                        max_column=max_column,
                        path=path + [key],
                    )
        elif _check_model_class(annotation, typing.Collection):
            for s in filter(lambda c: issubclass(c, BaseModel), typing.get_args(annotation)):
                _annotate_model(
                    s, comment_map[key], level=level + 1, max_column=max_column, path=path + [key]
                )


def _is_union_type(typ_) -> bool:
    return typing.get_origin(typ_) in (typing.Union, types.UnionType)


def _placeholder_model_instance(model: type[BaseModel], placeholder=None):
    """
    Constructs a model instance where the values of required fields are replaced by a placeholder
    """
    placeholders = {}
    # recurse into (optional) fields which are themselves pydantic models
    for name, field in model.model_fields.items():
        annotation = field.annotation

        # optional fields, i.e. `Union[Model, None]` or `Model | None`
        if _is_union_type(annotation):
            if len(args := typing.get_args(annotation)) == 2 and types.NoneType in args:
                if field.default in (PydanticUndefined, None):
                    sub_model = next(filter(lambda s: s is not types.NoneType, args))
                    if _check_model_class(sub_model):
                        placeholders[name] = _placeholder_model_instance(sub_model, placeholder)

        # required fields, i.e. `Model`
        if _check_model_class(annotation):
            placeholders[name] = _placeholder_model_instance(annotation)

    # replace values of undefined required fields with `placeholder`
    required_field_placeholders = {
        name: placeholder
        for name, field in model.model_fields.items()
        if field.is_required() and field.default is PydanticUndefined
    }
    required_field_placeholders.update(placeholders)

    # construct a model instance with invalid data, skipping validation!
    # this is only to be used to generate example configuration files
    invalid_model_instance = model.model_construct(_fields_set=None, **required_field_placeholders)
    return invalid_model_instance


def _yaml_instance():
    yaml = YAML(typ="rt")
    yaml.indent(mapping=INDENTATION, sequence=INDENTATION * 2, offset=INDENTATION)
    return yaml


def _load_yaml(yaml_str: str) -> ruamel.yaml.CommentedMap:
    return _yaml_instance().load(yaml_str)


def _dump_yaml(comment_map: ruamel.yaml.CommentedMap) -> str:
    yaml = _yaml_instance()

    with StringIO() as out:
        yaml.dump(comment_map, stream=out)
        return out.getvalue()


def _dump_commented_yaml(model: type[BaseModel], comment_optional: bool = True) -> str:
    invalid_model_instance = _placeholder_model_instance(model)

    cfg, max_column = _model_to_commented_yaml(invalid_model_instance)
    max_column = max(50, max_column)
    _annotate_model(model, cfg, max_column=max_column)
    key_paths = _optional_key_paths(model, cfg)
    cfg_yaml = _dump_yaml(cfg)
    if comment_optional:
        return _comment_key_paths_naive(cfg_yaml, key_paths)
    else:
        return cfg_yaml


def _model_to_commented_yaml(model_instance: BaseModel, **kwargs):
    yaml = _yaml_instance()
    with StringIO() as s:
        yaml.dump(
            json.loads(model_instance.model_dump_json(warnings=False, **kwargs)),
            stream=s,
        )
        s.flush()
        yaml_config_string = s.getvalue()
        max_column = max(map(len, yaml_config_string.splitlines())) + 2
        cfg = yaml.load(stream=yaml_config_string)
    return cfg, max_column


def _comment_key_paths_naive(
    yaml_str: str, key_paths: list[list[str]], comment_prefix: str = "#"
) -> str:
    comment_lines: set[int] = set()

    for key_path in key_paths:
        same_block = False
        key_line_indent = 0
        key = key_path[0]
        for i, line in enumerate(yaml_str.splitlines()):
            line_indent = len(line) - len(line.lstrip())
            key_match = line.lstrip().startswith(key + ":")
            if key_match:
                key_path.pop(0)
                if key_path:
                    key = key_path[0]
                    continue
                else:
                    key_line_indent = line_indent
                    same_block = True
                    comment_lines.add(i)
                    continue
            if same_block:
                if line_indent > key_line_indent or len(line) == 0:
                    comment_lines.add(i)
                else:
                    break
            else:
                if not key_path:
                    break

    return "\n".join(
        (comment_prefix + line) if i in comment_lines else line
        for i, line in enumerate(yaml_str.splitlines())
    )


def _optional_key_paths(
    config_model: type[BaseModel],
    comment_map: ruamel.yaml.CommentedMap,
    path: list[str] = [],
):
    optional_keys = []
    for key, field in config_model.model_fields.items():
        path_ = path + [key]
        annotation = field.annotation
        is_union = _is_union_type(annotation)
        allows_none = types.NoneType in typing_extensions.get_args(annotation)
        is_optional = (is_union and allows_none) or not field.is_required()
        if is_optional:
            optional_keys.append(path_)

        if _check_model_class(annotation):
            optional_keys.extend(_optional_key_paths(annotation, comment_map[key], path_))

    return optional_keys


class ToggleModel(SnappyModel):
    enabled: bool = False

    def __init__(self, enabled: bool = False, **kwargs):
        super().__init__(enabled=enabled, **kwargs)


__all__ = [
    "RelationshipDefinition",
    "SnappyModel",
    "SnappyStepModel",
    "ToggleModel",
]
