# -*- coding: utf-8 -*-
"""Tests for the ASCAT params of ``somatic_purity_ploidy_estimate``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.somatic_purity_ploidy_estimate import AscatStepPart


@pytest.mark.parametrize("action", ["cnv_tumor", "cnv_normal"])
def test_ascat_build_cnv_takes_no_params(action):
    # run_ascat names the sample column after the tumor library
    assert not hasattr(AscatStepPart, f"_get_params_{action}")


def test_ascat_run_names_the_tumor_library():
    part = object.__new__(AscatStepPart)

    params = part._get_params_run_ascat(SimpleNamespace(tumor_library="T1-DNA1-WGS1"))

    assert params == {"tumor_library": "T1-DNA1-WGS1"}
