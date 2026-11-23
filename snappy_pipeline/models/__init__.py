import enum
import os
import re
import textwrap
import types
import typing
from enum import Enum
from inspect import isclass
from io import StringIO
from typing import Annotated, ClassVar, Literal

import typing_extensions
from annotated_types import Predicate
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    ValidationInfo,
    model_validator,
)
from pydantic_core import PydanticUndefined
from ruamel.yaml import YAML
from snakemake.exceptions import UndefinedPathvarException

from snappy_pipeline.reads import compile_search_pattern


def _resolve_path(value: str, info: ValidationInfo) -> str:
    """Expand the ``pathvars`` in ``value`` and make it absolute against the project directory.

    ``load_project()`` gives both in the validation context, which also collects the paths for
    ``snappy check-files``; without that context the value is unchanged. Existence is not checked:
    Snakemake checks the inputs of the jobs it runs, and dry-runs work without the files.
    """
    context = info.context or {}
    if not value or "project_dir" not in context:
        return value
    try:
        expanded = context["pathvars"].apply(value)
    except UndefinedPathvarException as e:
        raise ValueError(str(e)) from e
    path = os.path.normpath(os.path.join(context["project_dir"], expanded))
    context["paths"].append(path)
    return path


ConfigPath = Annotated[str, AfterValidator(_resolve_path)]
"""A file or directory in the configuration: a task-local file or an ``external_data`` path."""

#: The scheme of an image URI, such as ``docker:`` or ``oras:``
_IMAGE_SCHEME = re.compile(r"[a-z][a-z0-9+.-]*:")


def _image(value: str, info: ValidationInfo) -> str:
    """Keep an image URI; resolve a local image file like a ``ConfigPath``."""
    return value if _IMAGE_SCHEME.match(value) else _resolve_path(value, info)


ContainerImage = Annotated[str, AfterValidator(_image)]
"""A container image: a URI such as ``docker://...``, or the path of a local ``.sif`` file."""


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


class SnappyModel(BaseModel):
    """
    Base class for all snappy models.
    By default, extra fields are forbidden, attribute docstrings are used for field descriptions,
    enum member values instead of names are used, and default values are validated (because
    validation can potentially modify the values of fields with default values)

    Path fields are ``ConfigPath``: ``load_project()`` expands their ``pathvars`` and resolves
    them against the project directory
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


class ToolModel(SnappyModel):
    """Base class of the config section of one tool

    Its ``extra_arguments`` reach the params of the step part's command (``tool_section``,
    ``extra_arguments_action`` of ``BaseStepPart``), and the wrapper puts them on the tool's
    command line.
    """

    extra_arguments: list[str] = []
    """Additional arguments for the tool's command line, such as ``["-C"]``"""


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

        # Find the index/proband library for any library of the same family
        index_lib:
          via: related_group
          target: "role == 'index'"
          column: index_lib

        # Find all tumor libraries in the same donor (one-to-many)
        donor_tumor_libs:
          via: donor_name
          target: "role == 'tumor'"
          column: donor_tumor_libs
          many: true
    """

    via: str | list[str]
    """Column name in the library DataFrame to use as the join key, or several names.
    The related row must have the *same* value in each of these columns as the
    source row, e.g. ``[donor_name, sample_name]`` for the libraries of one sample."""

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


class ConfigMode(enum.StrEnum):
    """Which fields the generated default configuration of a task shows
    (:func:`default_config_yaml_string`)."""

    COMMENTED = "commented"
    """Every field with its default; the optional ones commented out"""

    FULL = "full"
    """Every field with its default"""

    MINIMAL = "minimal"
    """The required fields only, with the ``depends_on`` fields that a tool requires; no
    descriptions, examples or options, only the ``REQUIRED`` tags"""


class Intermediate(enum.StrEnum):
    """What a task deletes once its readers are done (``delete_intermediates``)"""

    intermediates = "intermediates"
    """Per-sample intermediates that only later rules of the task read"""

    bam = "bam"
    """The BAM file of the task and its index, which another task processes further"""


class ReadBy:
    """``depends_on`` field metadata: who reads the field if it is set.

    A reader is a set of tools of the step (none: any tool) and an optional condition on a
    setting: ``when`` is a dotted path into the step's config, which is on if every toggle along
    it is enabled, every section on it is present and a bool at its end is true; ``value``
    compares the setting instead, ``contains`` tests membership in a list. A task must not set a
    field that none of its readers reads. A field without markers is read by every tool.
    A path that ends in an optional section (``igv``) is on when the section is set.
    """

    def __init__(self, *tools: str, when: str | None = None, value=None, contains=None):
        #: The tool values, such as ``"mutect2"``
        self.tools = frozenset(str(getattr(tool, "value", tool)) for tool in tools)
        self.when = when
        self.value = value
        self.contains = contains

    def condition(self) -> str:
        """Return the condition as text: ``mutect2.contamination``, ``intervals: targeted``."""
        if self.when is None:
            return ""
        if self.value is not None:
            return f"{self.when}: {getattr(self.value, 'value', self.value)}"
        if self.contains is not None:
            return f"{getattr(self.contains, 'value', self.contains)} in {self.when}"
        return self.when

    def __str__(self) -> str:
        tools = f"tool: {', '.join(sorted(self.tools))}" if self.tools else ""
        return " with ".join(part for part in (tools, self.condition()) if part)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self})"

    def applies(self, model: BaseModel) -> bool:
        """Whether this reader reads the field in the task config ``model``."""
        tool = getattr(model, "tool", None)
        if self.tools and str(getattr(tool, "value", tool)) not in self.tools:
            return False
        if self.when is None:
            return True
        setting = model
        for name in self.when.split("."):
            setting = getattr(setting, name)
            if setting is None or (isinstance(setting, ToggleModel) and not setting.enabled):
                return False
        if self.value is not None:
            return setting == self.value
        if self.contains is not None:
            return self.contains in setting
        return bool(setting)


class RequiredBy(ReadBy):
    """``depends_on`` field metadata: a reader that needs the field (see :class:`ReadBy`)."""


def reader_markers(field) -> list[ReadBy]:
    """Return the :class:`ReadBy` and :class:`RequiredBy` markers of a ``depends_on`` field;
    empty for a field that every tool reads."""
    return [marker for marker in field.metadata if isinstance(marker, ReadBy)]


class DependsOnError(ValueError):
    """Raised for ``depends_on`` fields that no reader of the task reads (``unread``) or that a
    reader needs (``missing``, by reader)."""

    def __init__(self, missing: dict[str, list[str]], unread: dict[str, list[ReadBy]]):
        self.missing = missing
        self.unread = unread
        problems = [
            f"{reader} needs {' and '.join(f'depends_on.{name}' for name in names)}"
            for reader, names in missing.items()
        ]
        problems += [
            f"depends_on.{name} is read only by {' or '.join(map(str, markers))}"
            for name, markers in unread.items()
        ]
        super().__init__("; ".join(problems))


# This exists to distinguish the task config models of workflow steps from other snappy
# specific models.
class SnappyStepModel(SnappyModel, object):
    """A base class for all workflow step configuration models.

    Every step gets ``library_selection``, ``group_by``, and ``relationships``
    fields.  Steps that do not use these fields leave them at their defaults
    (``None`` / empty).
    """

    library_selection: str | None = None
    """A pandas ``DataFrame.query()`` expression over the library table that chooses the
    libraries this task processes, such as ``data_set == 'batch2'`` or ``library_type == 'WES'``;
    a task never takes libraries that its step cannot process. Unset: the step's default (the
    tumor DNA libraries of cancer sheets and the DNA libraries of germline sheets, unless the
    step's docs page names another one). The columns are those of the library table (see the
    sample sheets section of the docs): ``library_name``, ``data_set``, ``kind``, ``role``,
    ``extraction_type``, ``library_type``, ``library_kit``, ``donor_name``, ``sex``,
    ``related_group`` and the custom fields of the sheet."""

    background_selection: str | None = None
    """For steps that use the libraries of background data sets (mapping, panels of normals, the
    SV and gCNV callers): a pandas ``DataFrame.query()`` expression that narrows them, such as
    ``data_set == 'cohort_1000g'``. Unset: all background libraries."""

    delete_intermediates: list[Intermediate] = []
    """What the task deletes once its readers are done: ``intermediates``, the per-sample
    intermediate files that only later rules of the task read (raw calls, coverage, panel inputs,
    gVCFs, ...); ``bam``, the BAM file of a task whose BAM another task processes further
    (BQSR), which then is no result file. Scatter fragments are always
    deleted."""

    #: Whether the step writes a BAM file that another task can process further, so that
    #: ``delete_intermediates`` may name ``bam``
    bam_intermediate: ClassVar[bool] = False

    job_group: str | None = None
    """Snakemake job group of the task's jobs (local rules stay local): dependent jobs of the tasks
    with the same group run as one cluster job, and ``--group-components <group>=<n>`` puts n
    independent ones into one job. Their requests add up, so the members should share a
    partition."""

    group_by: Literal["library", "related"] | None = None
    """``library``: one output per library; ``related``: one output per group of related
    libraries (``related_group``), the members of a pedigree on germline sample sheets and the
    libraries of a donor on cancer sample sheets, named after the group: the sample sheet's family
    id (the data set's ``pedigree_field``, else ``familyId``), else the index donor with ``FAM_``
    (germline sheets), ``PAT_`` (cancer sheets) or no prefix (generic sheets).
    Unset: the step's default (related for joint germline calling), else the granularity that
    upstream tasks declare, else ``library``."""

    relationships: dict[str, RelationshipDefinition] | None = None
    """Named relationships that add derived columns to the library table before
    ``library_selection`` is applied. Each key is the column name; the value gives ``via`` (the
    join columns, such as ``donor_name``), ``target`` (a query over the related rows), and
    optionally ``column`` and ``many``."""

    def tool_tag(self) -> str | None:
        """The tag that the ``TOOL`` placeholder of a ``depends_on`` requirement stands for: the
        task's tool."""
        return getattr(self, "tool", None)

    @model_validator(mode="after")
    def validate_bam_intermediate(self):
        if Intermediate.bam in self.delete_intermediates and not self.bam_intermediate:
            raise ValueError("delete_intermediates: the step writes no BAM file to delete")
        return self

    @model_validator(mode="after")
    def validate_tool_readers(self):
        """Reject ``depends_on`` fields that none of their readers reads, and require those that
        an active reader needs (:class:`ReadBy`, :class:`RequiredBy`)."""
        depends_on = getattr(self, "depends_on", None)
        if depends_on is None:
            return self
        tool = getattr(self, "tool", None)
        missing, unread = {}, {}
        for name, field in type(depends_on).model_fields.items():
            markers = reader_markers(field)
            if not markers:
                continue
            active = [marker for marker in markers if marker.applies(self)]
            if getattr(depends_on, name):
                if not active:
                    unread[name] = markers
            elif needs := [marker for marker in active if isinstance(marker, RequiredBy)]:
                marker = needs[0]
                reader = f"tool={getattr(tool, 'value', tool)}" if marker.tools else ""
                reader = " with ".join(part for part in (reader, marker.condition()) if part)
                missing.setdefault(reader, []).append(name)
        if missing or unread:
            raise DependsOnError(missing, unread)
        return self

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
        cls,
        mode: ConfigMode | str = ConfigMode.COMMENTED,
        tool: str | None = None,
        answers: typing.Mapping[str, typing.Any] | None = None,
    ) -> str:
        """Return the annotated default ``config`` block of a task of this step
        (see :func:`default_config_yaml_string`)."""
        return default_config_yaml_string(cls, mode, tool, answers)


# Generation of the annotated default configuration of a step from its pydantic model, for
# ``snappy task add`` and the step pages of the docs: every field with its default, preceded by
# its description as a comment and followed by its tags (REQUIRED, the readers of a ``depends_on``
# field, options and examples) at a fixed column.

INDENTATION = 2  # Indentation used for YAML

#: Width of the generated YAML, comments included
LINE_WIDTH = 100

#: Column of the tags behind a key, unless a key line is longer
TAG_COLUMN = 50


def default_config_yaml_string(
    model: type[BaseModel],
    mode: ConfigMode | str = ConfigMode.COMMENTED,
    tool: str | None = None,
    answers: typing.Mapping[str, typing.Any] | None = None,
) -> str:
    """Return the annotated default configuration of ``model`` as YAML.

    A field's description precedes it as a comment, its tags follow it at :data:`TAG_COLUMN`.
    A blank line separates a described field or a section from the field before it. ``commented``
    comments out the optional fields and sections (a ``#`` right before the key). ``minimal``
    leaves them out, except for ``depends_on`` fields that a tool requires, and writes neither
    descriptions nor blank lines; of the tags, only the ``REQUIRED`` ones remain.

    With ``tool``, the block is that of a task with this tool: the sections of the other tools
    and the ``depends_on`` keys that only they read are left out, the tool's section is written,
    and a key the tool needs regardless of a setting is ``REQUIRED``. ``answers`` are values
    written in place of the defaults, nested like the block; an answered key is written without
    tags, an answered ``None`` as a key to fill in, tagged ``REQUIRED``.
    """
    answers = dict(answers or {})
    if tool is not None:
        answers.setdefault("tool", tool)
    writer = _ConfigWriter(ConfigMode(mode), tool)
    writer.write_model(model, _placeholder_data(model), answers, indent=0, commented=False)
    return writer.render()


def tool_values(model: type[BaseModel]) -> list[str]:
    """The tool values of a step model, from the enum or literal of its ``tool`` field;
    empty for a step without a tool."""
    field = model.model_fields.get("tool")
    if field is None:
        return []
    annotation = field.annotation
    if isinstance(annotation, type) and issubclass(annotation, enum.Enum):
        return [member.value for member in annotation]
    if typing.get_origin(annotation) is typing.Literal:
        return [str(value) for value in typing.get_args(annotation)]
    return []


class _ConfigWriter:
    """Collects the lines of the configuration with their tags, then aligns the tags."""

    def __init__(self, mode: ConfigMode, tool: str | None = None):
        self.mode = mode
        self.tool = tool
        #: (line, tags, width of the line without the comment char), so that the tags of the
        #: commented and the full configuration sit at the same column
        self.lines: list[tuple[str, str, int]] = []

    def write_model(
        self,
        model: type[BaseModel],
        data: dict,
        answers: typing.Mapping[str, typing.Any],
        indent: int,
        commented: bool,
    ):
        prefix = " " * indent
        documented = self.mode is not ConfigMode.MINIMAL
        sections = set(tool_values(model))
        first = True
        for name, field in model.model_fields.items():
            markers = reader_markers(field)
            if self.tool is not None:
                if name in sections and name != self.tool:
                    continue  # the section of another tool
                markers = [m for m in markers if not m.tools or self.tool in m.tools]
                if reader_markers(field) and not markers:
                    continue  # a key that only other tools read
            answered = name in answers
            needed = self.tool is not None and (
                name == self.tool
                or any(isinstance(m, RequiredBy) and m.when is None for m in markers)
            )
            required = not _is_optional(field) or needed or (answered and answers[name] is None)
            required_by = any(isinstance(m, RequiredBy) for m in markers)
            if not documented and not required and not required_by and not answered:
                continue
            comment_out = commented or (
                self.mode is ConfigMode.COMMENTED and not required and not answered
            )
            description = _description_lines(field, LINE_WIDTH - indent - 2) if documented else []
            sub_model = _sub_model(field.annotation)
            if documented and not first and (description or sub_model is not None):
                self.lines.append(("", "", 0))
            first = False
            self.lines += [(f"{prefix}# {line}".rstrip(), "", 0) for line in description]
            if (answered and answers[name] is not None) or name == self.tool:
                tags = ""  # an answered key, or the section of the tool, which is just written
            else:
                tags = _tags(field, documented, markers, required, self.tool)
            value = answers[name] if answered else data.get(name)
            if needed and value == "":  # a key the tool needs, without its empty default
                value = None
            if sub_model is not None:
                self._add(f"{prefix}{name}:", comment_out, tags)
                if not isinstance(value, dict):
                    value = _placeholder_data(sub_model)
                sub_answers = answers.get(name) if isinstance(answers.get(name), dict) else {}
                before = len(self.lines)
                self.write_model(sub_model, value, sub_answers, indent + INDENTATION, comment_out)
                if len(self.lines) == before:  # a section without lines is still a section
                    self.lines[-1] = (_commented(f"{prefix}{name}: {{}}", comment_out), tags, 0)
            elif answered and isinstance(value, dict) and None in value.values():
                # Answered keys to fill in, such as the files of an external_data task
                self._add(f"{prefix}{name}:", comment_out, tags)
                for sub_name, sub_value in value.items():
                    first_line, *rest = _yaml_lines(sub_name, sub_value)
                    sub_prefix = prefix + " " * INDENTATION
                    self._add(
                        sub_prefix + first_line,
                        comment_out,
                        "REQUIRED" if sub_value is None else "",
                    )
                    for line in rest:
                        self._add(sub_prefix + line, comment_out, "")
            else:
                first_line, *rest = _yaml_lines(name, value)
                self._add(prefix + first_line, comment_out, tags)
                for line in rest:
                    self._add(prefix + line, comment_out, "")

    def _add(self, line: str, commented: bool, tags: str) -> None:
        self.lines.append((_commented(line, commented), tags, len(line)))

    def render(self) -> str:
        if not self.lines:  # a minimal configuration of a step without required fields
            return ""
        column = max([TAG_COLUMN] + [width + 2 for _, tags, width in self.lines if tags])
        return (
            "\n".join(f"{text:<{column}}# {tags}" if tags else text for text, tags, _ in self.lines)
            + "\n"
        )


def _commented(line: str, commented: bool) -> str:
    """``line`` with a ``#`` right before its content, if ``commented``."""
    if not commented:
        return line
    content = line.lstrip()
    return line[: len(line) - len(content)] + "#" + content


def _yaml_lines(key: str, value: typing.Any) -> list[str]:
    """The YAML of ``key: value``: one line, or the key and the indented lines of a collection."""
    with StringIO() as out:
        _yaml_instance().dump({key: value}, out)
        return out.getvalue().rstrip("\n").splitlines()


def _description_lines(field, width: int) -> list[str]:
    """The field's description as comment text: paragraphs separated by an empty line, each
    wrapped to ``width``, list items with a hanging indent, literals in single backticks."""
    text = (field.description or "").strip().replace("``", "`")
    lines: list[str] = []
    for paragraph in re.split(r"\n\s*\n", text):
        if lines:
            lines.append("")
        for block in re.split(r"\n(?=- )", paragraph):
            joined = " ".join(line.strip() for line in block.splitlines())
            indent = "  " if joined.startswith("- ") else ""
            lines += textwrap.wrap(joined, width, subsequent_indent=indent)
    return lines


def _tags(
    field,
    documented: bool = True,
    markers: list[ReadBy] | None = None,
    required: bool | None = None,
    tool: str | None = None,
) -> str:
    """The tags behind a key: whether it is required and by which readers; if ``documented``,
    also its optional readers, examples and options. For a task with ``tool``, ``markers`` are
    the readers of that tool and ``required`` is settled, so a reader is named by its condition
    only, and one without a condition not at all."""
    if markers is None:
        markers = reader_markers(field)
    if required is None:
        required = not _is_optional(field)

    def reader(marker: ReadBy) -> str:
        return str(marker) if tool is None else marker.condition()

    tags = ["REQUIRED"] if required else []
    needs = [m for m in markers if isinstance(m, RequiredBy) and (tool is None or m.when)]
    if needs and not required:
        tags.append(f"REQUIRED for {'; '.join(reader(m) for m in needs)}")
    if not documented:
        return "; ".join(tags)
    reads = [m for m in markers if not isinstance(m, RequiredBy) and (tool is None or m.when)]
    if reads:
        tags.append(f"optional for {'; '.join(reader(m) for m in reads)}")
    if field.examples:
        tags.append(f"Examples: {', '.join(map(str, field.examples))}")
    extra = field.json_schema_extra if isinstance(field.json_schema_extra, dict) else {}
    options = extra.get("options", [])
    if _check_model_class(field.annotation, enum.Enum):
        options = enum_options(field.annotation)
    if options:
        tags.append("Options: " + ", ".join(_option(name, value) for name, value in options))
    return "; ".join(tags)


def _option(name: str, value: typing.Any) -> str:
    """``'value'``, with the member name if it differs: ``'local-smart' (LOCAL_SMART)``."""
    if name.upper() != str(value).upper().replace("-", "_"):
        return f"{value!r} ({name})"
    return repr(value)


def _is_optional(field) -> bool:
    """Whether the field has a default or accepts ``None``."""
    annotation = field.annotation
    allows_none = types.NoneType in typing_extensions.get_args(annotation)
    return (_is_union_type(annotation) and allows_none) or not field.is_required()


def _sub_model(annotation) -> type[BaseModel] | None:
    """The model class of a field typed ``Model`` or ``Model | None``, else ``None``."""
    if _check_model_class(annotation):
        return annotation
    if _is_union_type(annotation):
        args = typing.get_args(annotation)
        if len(args) == 2 and types.NoneType in args:
            other = next(arg for arg in args if arg is not types.NoneType)
            if _check_model_class(other):
                return other
    return None


def _check_model_class(annotation, clazz: type = BaseModel) -> bool:
    """Checks whether the given annotation is a class and is a subclass of `clazz`"""
    try:
        is_model_class = isclass(annotation) and issubclass(annotation, clazz)
    except TypeError:
        # for some reason, types.GenericAlias doesn't count as a class?
        is_model_class = hasattr(annotation, "model_fields")
    return is_model_class


def _is_union_type(typ_) -> bool:
    return typing.get_origin(typ_) in (typing.Union, types.UnionType)


def _placeholder_data(model: type[BaseModel]) -> dict:
    """The defaults of ``model`` as plain data, ``None`` for a required field without one."""
    return _placeholder_model_instance(model).model_dump(mode="json", warnings=False)


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


class LibraryFolders(SnappyModel):
    """Files in a folder per library below search paths: the reads of a data set, or the files of
    an ``external_data`` task in ``search_paths`` mode. ``BaseStep.read_groups()`` finds the reads
    of both the same way."""

    search_paths: list[ConfigPath] = []
    """Directories that contain the folder of each library, at any depth"""

    search_patterns: list[dict[str, str]] = []
    """Regular expressions for the paths below a library's folder, by key. Reads need a ``left``
    key (and ``right`` for paired-end data) and a ``(?P<readgroup>...)`` group, which pairs the
    mates of one read group."""

    mixed_se_pe: bool = False
    """For reads: accept read groups without a ``right`` file, so that single-end and paired-end
    data can be mixed"""

    def patterns_find_reads(self) -> bool:
        """Whether the patterns find reads, or one file per key."""
        return True

    @model_validator(mode="after")
    def check_search_patterns(self):
        if self.search_paths and not self.search_patterns:
            raise ValueError("search_paths need search_patterns")
        for pattern in self.search_patterns:
            compile_search_pattern(pattern, reads=self.patterns_find_reads())
        if self.mixed_se_pe and not self.patterns_find_reads():
            raise ValueError("mixed_se_pe is only for reads")
        return self


class ToggleModel(SnappyModel):
    enabled: bool = False

    def __init__(self, enabled: bool = False, **kwargs):
        super().__init__(enabled=enabled, **kwargs)


__all__ = [
    "ConfigPath",
    "ContainerImage",
    "LibraryFolders",
    "RelationshipDefinition",
    "SnappyModel",
    "SnappyStepModel",
    "ToggleModel",
]
