# -*- coding: utf-8 -*-
"""Tests for the CNVkit params of ``somatic_wgs_cnv_calling``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_wgs_cnv_calling import CnvkitSomaticWgsStepPart


def test_cnvkit_postprocess_gets_the_library_name():
    # The params do not depend on the config, so skip the full step setup.
    part = object.__new__(CnvkitSomaticWgsStepPart)

    params = part._get_params_postprocess(SimpleNamespace(library_name="T1-DNA1-WGS1"))

    assert params == {"library_name": "T1-DNA1-WGS1"}
