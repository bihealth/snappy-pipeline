# -*- coding: utf-8 -*-
"""Tests for the project configuration model ``ConfigModel`` and shared model validators."""

import re
from typing import Annotated

import pydantic
import pytest

from ruamel.yaml import YAML

from snappy_pipeline.models import (
    LINE_WIDTH,
    ConfigMode,
    ReadBy,
    RequiredBy,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    default_config_yaml_string,
)
from snappy_pipeline.workflow_model import ConfigModel, DataSet
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.gene_expression_quantification.model import (
    GeneExpressionQuantification,
)


def _minimal_config(**overrides):
    config = {
        "tasks": [{"step": "ngs_mapping", "name": "mapping", "config": {"tool": "bwa"}}],
        "data_sets": {},
    }
    config.update(overrides)
    return config


def test_config_model_accepts_minimal_config():
    config = ConfigModel(**_minimal_config())
    assert [t.name for t in config.tasks] == ["mapping"]
    # The step config stays a plain dict; the step's own model validates it later.
    assert config.tasks[0].config == {"tool": "bwa"}


@pytest.mark.parametrize("missing", ["tasks", "data_sets"])
def test_config_model_requires_top_level_sections(missing):
    config = _minimal_config()
    del config[missing]
    with pytest.raises(pydantic.ValidationError, match=missing):
        ConfigModel(**config)


def test_depends_on_belongs_into_the_task_config():
    task = {
        "step": "variant_calling",
        "name": "calling",
        "depends_on": {"alignments": "mapping"},
        "config": {},
    }
    with pytest.raises(pydantic.ValidationError, match="depends_on"):
        ConfigModel(**_minimal_config(tasks=[task]))


class _QcDependsOn(SnappyModel):
    reference: str
    reads: Annotated[str, RequiredBy("fastqc")] = ""
    dbsnp: Annotated[str, ReadBy("picard")] = ""
    kits: Annotated[
        str, RequiredBy("alfred"), RequiredBy("picard", when="picard.programs", contains="hs")
    ] = ""


class _Pileup(ToggleModel):
    min_depth: int = 1


class _Picard(SnappyModel):
    programs: list[str] = []
    pileup: _Pileup = _Pileup()


class _Qc(SnappyStepModel):
    depends_on: _QcDependsOn
    tool: str
    picard: _Picard | None = None


def test_tool_readers_require_and_reject_depends_on_fields():
    _Qc(tool="alfred", depends_on={"reference": "genome", "kits": "kits"})
    # picard reads dbsnp if it is set, and needs kits only for one program
    _Qc(tool="picard", depends_on={"reference": "genome"}, picard={})
    _Qc(tool="picard", depends_on={"reference": "genome", "dbsnp": "dbsnp"}, picard={})
    _Qc(tool="picard", depends_on={"reference": "genome", "kits": "k"}, picard={"programs": ["hs"]})

    with pytest.raises(pydantic.ValidationError) as error:
        _Qc(tool="fastqc", depends_on={"reference": "genome", "dbsnp": "dbsnp", "kits": "kits"})
    # All problems of the task in one message
    assert (
        "tool=fastqc needs depends_on.reads; depends_on.dbsnp is read only by tool: picard; "
        "depends_on.kits is read only by tool: alfred or tool: picard with hs in picard.programs"
    ) in str(error.value)

    with pytest.raises(pydantic.ValidationError, match="tool=picard with hs in picard.programs"):
        _Qc(tool="picard", depends_on={"reference": "genome"}, picard={"programs": ["hs"]})
    with pytest.raises(pydantic.ValidationError, match="depends_on.kits is read only by"):
        _Qc(tool="picard", depends_on={"reference": "genome", "kits": "kits"}, picard={})


class _SettingsDependsOn(SnappyModel):
    pileup_sites: Annotated[str, RequiredBy(when="picard.pileup")] = ""
    targets: Annotated[str, RequiredBy(when="intervals", value="targeted")] = ""
    any_flag: Annotated[str, ReadBy(when="flag"), ReadBy(when="picard.pileup")] = ""


class _Settings(SnappyStepModel):
    depends_on: _SettingsDependsOn
    intervals: str = "wgs"
    flag: bool = False
    picard: _Picard | None = None


@pytest.mark.parametrize(
    "settings, depends_on, error",
    [
        ({}, {}, None),
        # A missing section and a disabled toggle along the path are off
        ({}, {"pileup_sites": "x"}, "depends_on.pileup_sites is read only by picard.pileup"),
        ({"picard": {"pileup": {"enabled": False}}}, {"pileup_sites": "x"}, "read only by"),
        ({"picard": {"pileup": {"enabled": True}}}, {}, "picard.pileup needs depends_on.pileup"),
        ({"picard": {"pileup": {"enabled": True}}}, {"pileup_sites": "x", "any_flag": "y"}, None),
        # A value; a bool
        ({"intervals": "targeted"}, {}, "intervals: targeted needs depends_on.targets"),
        ({}, {"targets": "x"}, "depends_on.targets is read only by intervals: targeted"),
        ({"flag": True}, {"any_flag": "y"}, None),
        ({}, {"any_flag": "y"}, "read only by flag or picard.pileup"),
    ],
)
def test_tool_readers_with_settings(settings, depends_on, error):
    if error is None:
        _Settings(depends_on=depends_on, **settings)
    else:
        with pytest.raises(pydantic.ValidationError, match=error):
            _Settings(depends_on=depends_on, **settings)


def test_default_config_comments_sections_named_like_depends_on_keys():
    # gene_expression_quantification has depends_on.strandedness and the optional section of the
    # strandedness tool
    lines = GeneExpressionQuantification.default_config_yaml_string().splitlines()
    assert "#strandedness:" in lines
    assert "strandedness:" not in lines


def _uncommented(text: str) -> str:
    """The commented default configuration with the ``#`` before the keys removed; the tags
    behind the keys lose their alignment, so they are set off by one space."""
    lines = []
    for line in text.splitlines():
        content = line.lstrip()
        if content.startswith("#") and not content.startswith("# ") and content != "#":
            line = line[: len(line) - len(content)] + content[1:]
        lines.append(re.sub(r"\S\s+# ", lambda m: m.group()[0] + " # ", line))
    return "\n".join(lines) + "\n"


def _keys(data, path=()):
    if isinstance(data, dict):
        for key, value in data.items():
            yield from _keys(value, path + (key,))
    elif path:
        yield path


@pytest.mark.parametrize("step", sorted(WORKFLOW_REGISTRY))
def test_default_config_is_formatted_consistently(step):
    model = WORKFLOW_REGISTRY[step].config_model_class
    commented = default_config_yaml_string(model, ConfigMode.COMMENTED)
    full = default_config_yaml_string(model, ConfigMode.FULL)
    minimal = default_config_yaml_string(model, ConfigMode.MINIMAL)
    # The modes agree: uncommenting every key gives the full configuration, which is valid YAML
    # with every key of the minimal one
    assert _uncommented(commented) == _uncommented(full)
    yaml = YAML(typ="safe")
    full_keys = set(_keys(yaml.load(full)))
    assert set(_keys(yaml.load(minimal))) <= full_keys
    for line in commented.splitlines():
        content = line.lstrip()
        # Comments are the description of the key below (``# text``) or a commented-out key (``#key``)
        assert not content.startswith("##"), line
        assert not content.startswith("# ") or "#" not in line.split("# ", 1)[1][:1], line
        # Descriptions are wrapped to the width; tags (options, readers) and values (image
        # digests) may be longer
        assert len(line) <= LINE_WIDTH or not content.startswith("# "), line
    # A described key is set off by a blank line, a blank line never ends up as a bare comment
    assert "\n#\n" not in commented and "\n\n\n" not in commented
    # Literals are quoted with single backticks, no Sphinx roles or double backticks remain
    assert "``" not in commented and ":func:" not in commented and ":class:" not in commented
    # The minimal configuration has neither descriptions nor blank lines, and only REQUIRED tags
    for line in minimal.splitlines():
        assert line and not line.lstrip().startswith("#"), line
        assert "#" not in line or line.split("#", 1)[1].startswith(" REQUIRED"), line


def test_default_config_marks_required_fields_and_readers():
    lines = default_config_yaml_string(WORKFLOW_REGISTRY["ngs_mapping"].config_model_class)
    assert "\ndepends_on:" in lines
    assert "  #features: ''" in lines and "# REQUIRED for tool: star" in lines
    assert "\ntool:" in lines and "# REQUIRED; Options: 'bwa', 'bwa_mem2'" in lines
    # The minimal configuration keeps the fields that a tool requires, not the optional ones
    minimal = default_config_yaml_string(
        WORKFLOW_REGISTRY["ngs_mapping"].config_model_class, "minimal"
    )
    assert "features: ''" in minimal and "bwa:" not in minimal
    assert "# REQUIRED for tool: star" in minimal and "Options:" not in minimal


def test_data_set_defaults_to_r1_r2_patterns_but_not_to_search_paths():
    data_set = DataSet()
    assert data_set.search_paths == []
    assert data_set.search_patterns == [
        {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"}
    ]
