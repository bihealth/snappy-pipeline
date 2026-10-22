# -*- coding: utf-8 -*-
"""Tests for the CNVkit and PureCN params and the sequenza outputs of
``somatic_targeted_seq_cnv_calling``."""

from types import SimpleNamespace

import pytest
from snakemake.io import is_flagged

from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    CnvKitStepPart,
    PureCNStepPart,
    SequenzaStepPart,
)
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling.model import Cnvkit, PureCn


def _cnvkit_part(**cfg):
    # The params only read the CNVkit config section, so skip the full step setup.
    part = object.__new__(CnvKitStepPart)
    part.cfg = Cnvkit(**cfg)
    return part


@pytest.mark.parametrize("action", ["call", "plot"])
def test_gender_and_male_reference_are_passed_to_call_and_plot(action):
    params = _cnvkit_part(gender="male", male_reference=True)._cnvkit_params(action)

    assert params["gender"] == "male"
    assert params["male_reference"] is True


def test_cnvkit_postprocess_gets_the_tumor_library_and_the_ploidy():
    part = _cnvkit_part(call={"ploidy": 3})
    params = part.get_params("postprocess")(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"library_name": "T1-DNA1-WES1", "ploidy": 3}


def test_guessed_gender_and_female_reference_are_not_passed():
    params = _cnvkit_part()._cnvkit_params("call")

    assert "gender" not in params
    assert "male_reference" not in params


def _purecn_part():
    # The params only read the PureCN config section, so skip the full step setup.
    part = object.__new__(PureCNStepPart)
    part.config = SimpleNamespace(get={"purecn": PureCn()}.get)
    return part


def test_purecn_coverage_gets_only_the_config():
    params = _purecn_part()._get_params_coverage(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert set(params) == {"config"}
    assert params["config"]["seed"] == 1234567


def test_purecn_run_gets_the_config_and_the_tumor_library():
    params = _purecn_part()._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert set(params) == {"config", "library_name"}
    assert params["library_name"] == "T1-DNA1-WES1"


def test_sequenza_run_writes_its_report_directory():
    outputs = object.__new__(SequenzaStepPart)._get_output_files_run()
    assert outputs["report"] == "work/{tumor_library}/report"
    assert is_flagged(outputs["report"], "directory")
    assert "report" in SequenzaStepPart.result_keys["run"]
