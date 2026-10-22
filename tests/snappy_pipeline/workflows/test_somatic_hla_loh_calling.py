# -*- coding: utf-8 -*-
"""Tests for the LOHHLA params of ``somatic_hla_loh_calling``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_hla_loh_calling import LohhlaStepPart


def test_lohhla_passes_the_minimum_coverage():
    # The params only read the step config, so skip the full step setup.
    part = object.__new__(LohhlaStepPart)
    part.config = SimpleNamespace(min_coverage_filter=20)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WGS1"))

    assert params == {"tumor_library": "T1-DNA1-WGS1", "min_coverage_filter": 20}
