import enum
import json
import types
import typing
from inspect import isclass
from io import StringIO

import ruamel
import typing_extensions
from pydantic import BaseModel
from pydantic_core import PydanticUndefined
from ruamel.yaml import YAML


def is_union_type(typ_) -> bool:
    return typing.get_origin(typ_) in (typing.Union, types.UnionType)


def placeholder_model_instance(model: typing.Type[BaseModel], placeholder=None):
    """
    Constructs a model instance where the value of required fields are replaced by a placeholder
    """
    placeholders = {}
    # recurse into (optional) fields which are themselves pydantic models
    for name, field in model.model_fields.items():
        annotation = field.annotation

        # optional fields, i.e. `Union[Model, None]` or `Model | None`
        if is_union_type(annotation):
            if len(args := typing.get_args(annotation)) == 2 and types.NoneType in args:
                if field.default is PydanticUndefined:
                    sub_model = next(filter(lambda s: s is not types.NoneType, args))
                    if issubclass(sub_model, BaseModel):
                        placeholders[name] = placeholder_model_instance(sub_model, placeholder)

        # required fields, i.e. `Model`
        if isclass(annotation) and issubclass(annotation, BaseModel):
            placeholders[name] = placeholder_model_instance(annotation)

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


def dump_commented_yaml(model: typing.Type[BaseModel]) -> str:
    yaml = YAML(typ="rt")
    yaml.indent(mapping=2, offset=2)
    yaml.block_seq_indent = 2
    invalid_model_instance = placeholder_model_instance(model)

    with StringIO() as s:
        yaml.dump(json.loads(invalid_model_instance.model_dump_json()), stream=s)
        s.flush()
        yaml_config_string = s.getvalue()
        max_column = max(map(len, yaml_config_string.splitlines())) + 2
        cfg = yaml.load(stream=yaml_config_string)

    annotate_model(model, cfg, max_column=max_column)
    key_paths = optional_key_paths(model, cfg)
    commented_yaml = comment_key_paths_naive(dump_yaml(cfg), key_paths)
    cfg = yaml.load(commented_yaml)

    with StringIO() as out:
        yaml.dump(cfg, stream=out)
        return out.getvalue()


def dump_yaml(comment_map: ruamel.yaml.CommentedMap) -> str:
    yaml = YAML(typ="rt")
    yaml.indent(mapping=2, offset=2)
    yaml.block_seq_indent = 2

    with StringIO() as out:
        yaml.dump(comment_map, stream=out)
        return out.getvalue()


def load_yaml(yaml_str: str) -> ruamel.yaml.CommentedMap:
    yaml = YAML(typ="rt")
    yaml.indent(mapping=2, offset=2)
    yaml.block_seq_indent = 2
    return yaml.load(yaml_str)


def comment_key_paths_naive(
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


def annotate_model(
    config_model: typing.Type[BaseModel],
    comment_map: ruamel.yaml.CommentedMap,
    max_column: int = 80,
    level: int = 0,
    indent: int = 2,
    path: list[str] = [],
):
    """
    Annotates a given pydantic model with comments
    based on whether a field is required or not, as well as listing options for enum valued fields
    """
    for key, field in config_model.model_fields.items():
        annotation = field.annotation
        is_union = is_union_type(annotation)
        allows_none = types.NoneType in typing_extensions.get_args(annotation)
        is_optional = (is_union and allows_none) or not field.is_required()

        if is_optional:
            tags = []
        else:
            tags = ["REQUIRED"]

        comment = tags

        if field.examples:
            comment.append(f"Examples: {', '.join(map(str, field.examples))}")

        options = []
        if field.json_schema_extra and "options" in field.json_schema_extra:
            options = field.json_schema_extra["options"]
        if isclass(annotation) and issubclass(annotation, enum.Enum):
            options = [(e.name, e.value) for e in annotation]

        if options:

            def option_str(option: tuple[typing.Any, typing.Any]) -> str:
                name, value = option
                if name != value:
                    return f"{value} ({name})"
                else:
                    return f"{value}"

            comment.append(f"Options: " f"{', '.join(map(option_str, options))}")

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

        if isclass(annotation) and issubclass(annotation, BaseModel):
            annotate_model(
                annotation,
                comment_map[key],
                level=level + 1,
                max_column=max_column,
                path=path + [key],
            )
        elif is_union_type(annotation):
            if len(args := typing.get_args(annotation)) == 2 and types.NoneType in args:
                sub_model = next(filter(lambda s: s is not types.NoneType, args))
                if isclass(sub_model) and issubclass(sub_model, BaseModel):
                    annotate_model(
                        sub_model,
                        comment_map[key],
                        level=level + 1,
                        max_column=max_column,
                        path=path + [key],
                    )
        elif isclass(annotation) and issubclass(annotation, typing.Collection):
            for s in filter(lambda c: issubclass(c, BaseModel), typing.get_args(annotation)):
                annotate_model(
                    s, comment_map[key], level=level + 1, max_column=max_column, path=path + [key]
                )


def optional_key_paths(
    config_model: typing.Type[BaseModel],
    comment_map: ruamel.yaml.CommentedMap,
    path: list[str] = [],
):
    optional_keys = []
    for key, field in config_model.model_fields.items():
        path_ = path + [key]
        annotation = field.annotation
        is_union = is_union_type(annotation)
        allows_none = types.NoneType in typing_extensions.get_args(annotation)
        is_optional = (is_union and allows_none) or not field.is_required()
        if is_optional:
            optional_keys.append(path_)

        if isclass(annotation) and issubclass(annotation, BaseModel):
            optional_keys.extend(optional_key_paths(annotation, comment_map[key], path_))

    return optional_keys
