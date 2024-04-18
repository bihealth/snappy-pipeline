import json
import typing
from copy import deepcopy
from io import StringIO

import pydantic
import ruamel
import typing_extensions
from ruamel.yaml import YAML


def placeholder_model_instance(model: typing.Type[pydantic.BaseModel], placeholder=None):
    """
    Constructs a model instance where the value of required fields are replaced by a placeholder
    """
    print("Doing stuff")
    placeholders = {}
    for name, field in model.model_fields.items():
        print(name)
        if issubclass(type(field), pydantic.BaseModel) or isinstance(field, pydantic.BaseModel):
            print(field.__class_)
            placeholders[name] = placeholder_model_instance(field.__class_)

    required_field_placeholders = {name: placeholder for name, field in
                                   model.model_fields.items() if
                                   field.is_required()}
    required_field_placeholders.update(placeholders)
    invalid_model_instance = model.model_construct(_fields_set=None,
                                                   **required_field_placeholders)
    return invalid_model_instance


def dump_commented_yaml(model: typing.Type[pydantic.BaseModel]) -> str:
    yaml = YAML()
    yaml.indent(mapping=2, offset=2)
    invalid_model_instance = placeholder_model_instance(model)

    with StringIO() as s:
        yaml.dump(json.loads(invalid_model_instance.model_dump_json()), stream=s)
        s.flush()
        yaml_config_string = s.getvalue()
        max_column = max(map(len, yaml_config_string.splitlines())) + 2
        cfg = yaml.load(stream=yaml_config_string)

    annotate_model(model, cfg, max_column=max_column)

    with StringIO() as out:
        yaml.dump(cfg, stream=out)
        return out.getvalue()


def annotate_model(config_model: typing.Type[pydantic.BaseModel],
                   comment_map: ruamel.yaml.CommentedMap,
                   max_column: int = 80,
                   level: int = 0,
                   indent: int = 2,
                   path: list[str] = []):
    for name, field in config_model.model_fields.items():
        print(name, "\t", field, field.description)
        key = name
        is_union = typing.Union == typing_extensions.get_origin(field.annotation)
        allows_none = type(None) in typing_extensions.get_args(field.annotation)

        if (is_union and allows_none) or not field.is_required():
            tags = []
        else:
            tags = ["REQUIRED"]

        comment = tags

        if field.examples:
            comment.append(f"Examples: {', '.join(map(str, field.examples))}")

        if options := (field.json_schema_extra or {}).get("options", None):
            comment.append(f"Options: {', '.join(map(str, options))}")

        comment = '; '.join(comment)

        if comment_map is not None:
            if comment:
                comment_map.yaml_add_eol_comment(comment, key, column=max_column)

            if field.description:
                (comment_map.
                 yaml_set_comment_before_after_key(key,
                                                   indent=indent * level,
                                                   before="\n" + field.description,
                                                   after=None, ))
        else:
            continue

        if (submodel := field.default) is not None:
            if isinstance(submodel, typing.Collection):
                for s in submodel:
                    if isinstance(s, pydantic.BaseModel):
                        print("annotating multiple", s)
                        annotate_model(s, comment_map[key], level=level + 1, max_column=max_column,
                                       path=path + [name])
            elif isinstance(submodel, pydantic.BaseModel):
                print("annotating submodel", type(submodel), submodel)
                annotate_model(submodel, comment_map[key], level=level + 1, max_column=max_column,
                               path=path + [name])
