# -*- coding: utf-8 -*-
"""Tests for the project configuration model ``ConfigModel`` and shared model validators."""

from typing import Annotated

import pydantic
import pytest

from snappy_pipeline.models import ReadBy, RequiredBy, SnappyModel, SnappyStepModel
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
    kits: Annotated[str, RequiredBy("alfred"), ReadBy("picard")] = ""


class _Qc(SnappyStepModel):
    depends_on: _QcDependsOn
    tool: str


def test_tool_readers_require_and_reject_depends_on_fields():
    _Qc(tool="alfred", depends_on={"reference": "genome", "kits": "kits"})
    # picard reads kits and dbsnp if they are set
    _Qc(tool="picard", depends_on={"reference": "genome"})
    _Qc(tool="picard", depends_on={"reference": "genome", "dbsnp": "dbsnp", "kits": "kits"})

    with pytest.raises(pydantic.ValidationError) as error:
        _Qc(tool="fastqc", depends_on={"reference": "genome", "dbsnp": "dbsnp", "kits": "kits"})
    # All problems of the task in one message
    assert (
        "tool=fastqc needs depends_on.reads; tool=fastqc does not read depends_on.dbsnp (read "
        "only by picard) and depends_on.kits (read only by alfred, picard)"
    ) in str(error.value)


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
