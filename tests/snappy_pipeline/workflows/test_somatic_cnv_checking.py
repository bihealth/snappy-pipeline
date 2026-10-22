# -*- coding: utf-8 -*-
"""Tests for the params of ``somatic_cnv_checking``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_cnv_checking import SomaticCnvCheckingCnvStepPart


def test_cnv_names_the_tumor_library():
    # The params do not depend on the config, so skip the full step setup.
    part = object.__new__(SomaticCnvCheckingCnvStepPart)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"tumor_library": "T1-DNA1-WES1"}
