# -*- coding: utf-8 -*-
"""Tests for the CNVkit and PureCN params of ``panel_of_normals``."""

from types import SimpleNamespace

from snappy_pipeline.models.cnvkit import PanelOfNormals as CnvKit
from snappy_pipeline.workflows.panel_of_normals import CnvkitStepPart, PureCnStepPart
from snappy_pipeline.workflows.panel_of_normals.model import PureCn


def test_purecn_coverage_gets_only_the_config_without_files():
    # The params only read the PureCN config section, so skip the full step setup.
    part = object.__new__(PureCnStepPart)
    part.config = SimpleNamespace(get={"purecn": PureCn(mappability="mappability.bw")}.get)

    params = part._get_params_coverage(SimpleNamespace(library_name="N1-DNA1-WES1"))

    assert set(params) == {"config"}
    assert "mappability" not in params["config"]


def test_cnvkit_target_passes_the_default_average_size():
    # The params only read the cnvkit config section, so skip the full step setup.
    part = object.__new__(CnvkitStepPart)
    part.name = "cnvkit"
    part.is_wgs = True
    part.config = SimpleNamespace(get=lambda name: CnvKit(target={"default_avg_size": 3000}))

    params = part._get_params_target(SimpleNamespace())

    assert params["default_avg_size"] == 3000
    assert CnvKit().target.default_avg_size == 5000


def test_purecn_container_is_a_config_field():
    assert PureCn().container == "docker://markusriester/purecn"
    assert PureCn(container="docker://markusriester/purecn:2.8.0").container.endswith(":2.8.0")
