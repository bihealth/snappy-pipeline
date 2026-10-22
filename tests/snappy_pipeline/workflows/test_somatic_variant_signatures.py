# -*- coding: utf-8 -*-
"""Tests for the VCF tabulation params of ``somatic_variant_signatures``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_variant_signatures import TabulateVariantsStepPart


def test_tabulate_vcf_names_the_tumor_also_without_matched_normal():
    # The params do not depend on the sample sheet, so skip the full step setup.
    part = object.__new__(TabulateVariantsStepPart)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"tumor_library": "T1-DNA1-WES1"}
