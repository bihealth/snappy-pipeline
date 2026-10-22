# -*- coding: utf-8 -*-
"""Tests for the params of ``somatic_cnv_checking``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_cnv_checking import (
    SomaticCnvCheckingCnvStepPart,
    SomaticCnvCheckingPileupStepPart,
    SomaticCnvCheckingReportStepPart,
)
from snappy_pipeline.workflows.somatic_cnv_checking.model import SomaticCnvChecking


def test_cnv_names_the_tumor_library():
    # The params do not depend on the config, so skip the full step setup.
    part = object.__new__(SomaticCnvCheckingCnvStepPart)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"tumor_library": "T1-DNA1-WES1"}


def _pileup_part():
    # The params only read the step config, so skip the full step setup.
    part = object.__new__(SomaticCnvCheckingPileupStepPart)
    part.config = SomaticCnvChecking.model_construct()
    return part


def test_normal_pileup_selects_heterozygous_variants_by_baf_and_depth():
    params = _pileup_part()._get_params_normal(SimpleNamespace(tumor_library="N1-DNA1-WES1"))

    assert params == {"min_baf": 0.4, "min_depth": 20, "max_depth": 10000}


def test_tumor_pileup_takes_the_sites_of_the_normal_without_filters():
    params = _pileup_part()._get_params_tumor(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"max_depth": 10000}


def test_report_passes_the_amplification_threshold():
    part = object.__new__(SomaticCnvCheckingReportStepPart)
    part.config = SomaticCnvChecking.model_construct(amplification=9)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"tumor_library": "T1-DNA1-WES1", "amplification": 9}
    assert SomaticCnvChecking.model_fields["amplification"].default == 7
