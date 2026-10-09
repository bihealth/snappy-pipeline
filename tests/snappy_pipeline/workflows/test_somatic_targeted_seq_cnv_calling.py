# -*- coding: utf-8 -*-
"""Tests for the CNVkit params of ``somatic_targeted_seq_cnv_calling``."""

import pytest

from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import CnvKitStepPart
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling.model import Cnvkit


def _cnvkit_part(**cfg):
    # The params only read the CNVkit config section, so skip the full step setup.
    part = object.__new__(CnvKitStepPart)
    part.cfg = Cnvkit(path_target="targets.bed", path_antitarget="antitargets.bed", **cfg)
    return part


@pytest.mark.parametrize("action", ["call", "plot"])
def test_gender_and_male_reference_are_passed_to_call_and_plot(action):
    params = _cnvkit_part(gender="male", male_reference=True)._cnvkit_params(action)

    assert params["gender"] == "male"
    assert params["male_reference"] is True


def test_guessed_gender_and_female_reference_are_not_passed():
    params = _cnvkit_part()._cnvkit_params("call")

    assert "gender" not in params
    assert "male_reference" not in params
