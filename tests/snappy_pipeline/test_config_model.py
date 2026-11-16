# -*- coding: utf-8 -*-
"""Tests for the project configuration model ``ConfigModel`` and shared model validators."""

from typing import Annotated

import pydantic
import pytest

from snappy_pipeline.models import (
    ReadBy,
    RequiredBy,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
)
from snappy_pipeline.workflow_model import ConfigModel, DataSet
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


def test_data_set_defaults_to_r1_r2_patterns_but_not_to_search_paths():
    data_set = DataSet()
    assert data_set.search_paths == []
    assert data_set.search_patterns == [
        {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"}
    ]
