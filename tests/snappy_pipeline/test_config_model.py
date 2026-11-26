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
from snappy_pipeline.workflow_model import ConfigModel
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.gene_expression_quantification.model import (
    GeneExpressionQuantification,
)


def _minimal_config(**overrides):
    config = {
        "project": {"name": "minimal"},
        "tasks": [{"step": "ngs_mapping", "name": "mapping", "config": {"tool": "bwa"}}],
        "sample_sheets": [],
    }
    config.update(overrides)
    return config


def test_config_model_accepts_minimal_config():
    config = ConfigModel(**_minimal_config())
    assert [t.name for t in config.tasks] == ["mapping"]
    # The step config stays a plain dict; the step's own model validates it later.
    assert config.tasks[0].config == {"tool": "bwa"}


@pytest.mark.parametrize("missing", ["project", "tasks"])
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


def test_default_config_of_a_tool_leaves_out_what_the_tool_does_not_read():
    model = WORKFLOW_REGISTRY["ngs_mapping"].config_model_class

    bwa = default_config_yaml_string(model, "minimal", tool="bwa")

    # The tool is set, its section written as such, the others and their keys left out
    assert "tool: bwa\n" in bwa and "bwa: {}\n" in bwa and "REQUIRED; Options" not in bwa
    assert "star" not in bwa and "features" not in bwa and "exons" not in bwa
    assert re.search(r"\n  index:\s+# REQUIRED\n", bwa)
    # A key the tool needs is REQUIRED, not "for" the tool
    star = default_config_yaml_string(model, "minimal", tool="star")
    assert re.search(r"\n  features:\s+# REQUIRED\n", star) and "features: ''" not in star
    assert "REQUIRED for" not in star and "\nstar: {}" in star
    commented = default_config_yaml_string(model, "commented", tool="bwa")
    assert "\nbwa:\n" in commented and "#bwa:" not in commented and "#star:" not in commented
    assert "#bwa_mem2:" not in commented and "REQUIRED for tool" not in commented
    # A key a setting needs keeps its condition, a key the tool needs is REQUIRED
    report = WORKFLOW_REGISTRY["variant_report"].config_model_class
    datavzrd = default_config_yaml_string(report, "minimal", tool="datavzrd")
    assert "# REQUIRED for igv" in datavzrd and "REQUIRED for tool" not in datavzrd
    filtration = WORKFLOW_REGISTRY["variant_filtration"].config_model_class
    dkfz = default_config_yaml_string(filtration, "minimal", tool="dkfz")
    assert re.search(r"\n  alignments:\s+# REQUIRED\n", dkfz) and "REQUIRED for" not in dkfz


def test_default_config_writes_the_answers_without_tags():
    external = WORKFLOW_REGISTRY["external_data"].config_model_class
    answers = {"produces": {"type": "raw", "tags": ["dna", "reference"]}, "files": {"fasta": None}}

    text = default_config_yaml_string(external, "minimal", answers=answers)

    # An answered key has no tag; an answered None is a key to fill in
    assert text.startswith("produces:\n  type: raw\n  tags:\n    - dna\n    - reference\n")
    assert re.search(r"\nfiles:\n  fasta:\s+# REQUIRED\n", text)
    mapping = WORKFLOW_REGISTRY["ngs_mapping"].config_model_class
    text = default_config_yaml_string(
        mapping, "commented", tool="bwa", answers={"depends_on": {"reads": "reads"}}
    )
    assert "\n  reads: reads\n" in text and re.search(r"\n  index:\s+# REQUIRED", text)


def test_sample_sheets_default_to_the_sheet_next_to_the_config():
    config = {"project": {"name": "x"}, "tasks": []}
    assert ConfigModel.model_validate(config).sample_sheets == ["samplesheet.tsv"]
    assert ConfigModel.model_validate(config | {"sample_sheets": []}).sample_sheets == []
