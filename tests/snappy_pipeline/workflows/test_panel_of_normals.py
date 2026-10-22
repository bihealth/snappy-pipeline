# -*- coding: utf-8 -*-
"""Tests for the PureCN params of ``panel_of_normals``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.panel_of_normals import PureCnStepPart
from snappy_pipeline.workflows.panel_of_normals.model import PureCn


def test_purecn_coverage_gets_only_the_config_without_files():
    # The params only read the PureCN config section, so skip the full step setup.
    part = object.__new__(PureCnStepPart)
    part.config = SimpleNamespace(get={"purecn": PureCn(mappability="mappability.bw")}.get)

    params = part._get_params_coverage(SimpleNamespace(library_name="N1-DNA1-WES1"))

    assert set(params) == {"config"}
    assert "mappability" not in params["config"]
