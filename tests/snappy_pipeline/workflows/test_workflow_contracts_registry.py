# -*- coding: utf-8 -*-
"""Contract tests for workflow registry definitions."""

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract import BaseStep


def test_registry_has_unique_keys_and_classes():
    assert len(WORKFLOW_REGISTRY) == len(set(WORKFLOW_REGISTRY.keys()))
    assert len(WORKFLOW_REGISTRY.values()) == len(set(WORKFLOW_REGISTRY.values()))


@pytest.mark.parametrize("step_name, workflow_cls", sorted(WORKFLOW_REGISTRY.items()))
def test_registry_entry_contract(step_name, workflow_cls):
    assert issubclass(workflow_cls, BaseStep)
    assert workflow_cls.name == step_name
    assert callable(getattr(workflow_cls, "default_config_yaml", None))


@pytest.mark.parametrize("step_name, workflow_cls", sorted(WORKFLOW_REGISTRY.items()))
def test_default_config_yaml_has_valid_top_level_shape(step_name, workflow_cls):
    _ = step_name
    yaml = ruamel_yaml.YAML()
    config_text = workflow_cls.default_config_yaml()
    assert isinstance(config_text, str)
    if not config_text.strip():
        # link_in is a pure config-carrier helper and currently ships an empty default snippet.
        assert step_name == "link_in"
        return

    config = yaml.load(config_text)
    assert isinstance(config, dict)
    assert "tasks" in config or "step_config" in config
