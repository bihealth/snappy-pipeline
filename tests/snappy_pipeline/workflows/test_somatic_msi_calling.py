# -*- coding: utf-8 -*-
"""Tests for the MANTIS params of ``somatic_msi_calling``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_msi_calling import Mantis2StepPart
from snappy_pipeline.workflows.somatic_msi_calling.model import SomaticMsiCalling

THRESHOLDS = {
    "min_read_length": 35,
    "min_read_quality": 20.0,
    "min_locus_quality": 25.0,
    "min_locus_coverage": 20,
    "min_repeat_reads": 1,
}


def test_mantis_thresholds_default_to_the_recommended_values():
    defaults = {name: SomaticMsiCalling.model_fields[name].default for name in THRESHOLDS}
    assert defaults == THRESHOLDS


def test_mantis_passes_its_thresholds():
    # The params only read the config, so skip the full step setup.
    part = object.__new__(Mantis2StepPart)
    part.config = SimpleNamespace(**{**THRESHOLDS, "min_locus_coverage": 30})

    params = part._get_params_run(SimpleNamespace(tumor_library="T1"))

    assert params == {**THRESHOLDS, "min_locus_coverage": 30}
