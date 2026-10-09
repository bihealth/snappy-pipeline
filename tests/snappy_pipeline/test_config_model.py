# -*- coding: utf-8 -*-
"""Tests for the project configuration model ``ConfigModel``."""

import pydantic
import pytest

from snappy_pipeline.workflow_model import ConfigModel


def _minimal_config(**overrides):
    config = {
        "static_data_config": {"reference": {"path": "/refs/genome.fa"}},
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


@pytest.mark.parametrize("missing", ["static_data_config", "tasks", "data_sets"])
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
