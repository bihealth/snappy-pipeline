# -*- coding: utf-8 -*-
"""Tests for the project configuration model ``ConfigModel`` and shared model validators."""

from types import SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.models import validators
from snappy_pipeline.models.cnvkit import Segment, Segmetrics
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


def test_require_tool_dependencies_names_missing_fields():
    required = {"picard": ("alignments", "reads")}
    model = SimpleNamespace(tool="picard", depends_on=SimpleNamespace(alignments="", reads=""))
    with pytest.raises(
        ValueError, match="tool=picard needs depends_on.alignments and depends_on.reads"
    ):
        validators.require_tool_dependencies(model, required)

    model.depends_on.alignments = model.depends_on.reads = "mapping"
    validators.require_tool_dependencies(model, required)
    # Tools without an entry need nothing.
    validators.require_tool_dependencies(SimpleNamespace(tool="fastqc"), required)


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


@pytest.mark.parametrize("threshold", [0, 2, 5])
def test_cnvkit_hmm_takes_a_smoothing_window(threshold):
    assert Segment(method="hmm-tumor", threshold=threshold).threshold == threshold


@pytest.mark.parametrize("threshold", [0.000001, 1, 2.5])
def test_cnvkit_hmm_rejects_a_threshold_that_is_no_window(threshold):
    with pytest.raises(pydantic.ValidationError, match="smoothing window"):
        Segment(method="hmm", threshold=threshold)


def test_cnvkit_smooth_bootstrap_is_a_bin_count_not_a_flag():
    assert Segmetrics(smooth_bootstrap=10).smooth_bootstrap == 10
    with pytest.raises(pydantic.ValidationError):
        Segmetrics(smooth_bootstrap=True)
