"""Tests for the FDR control of varlociraptor calls in ``variant_filtration``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_filtration import (
    VariantFiltrationWorkflow,
    VarlociraptorStepPart,
)
from snappy_pipeline.workflows.variant_filtration.model import VariantFiltration, VarlociraptorFdr


def test_the_fdr_params_are_the_settings_without_a_filter_name():
    part = object.__new__(VarlociraptorStepPart)
    part.config = SimpleNamespace(
        tool="varlociraptor", varlociraptor=VarlociraptorFdr(events=["present"], fdr=0.1)
    )
    part.parent = SimpleNamespace(task_name="fdr")

    params = part._get_params_run(SimpleNamespace())

    assert params == {
        "extra_arguments": [],
        "events": ["present"],
        "fdr": 0.1,
        "fdr_mode": "local-smart",
        "mode": "filter",
        "variant_types": ["SNV", "MNV", "INS", "DEL", "BND", "INV", "DUP", "REP"],
    }


def test_the_tag_mode_passes_the_task_name_as_filter_id():
    part = object.__new__(VarlociraptorStepPart)
    part.config = SimpleNamespace(
        tool="varlociraptor", varlociraptor=VarlociraptorFdr(events=["present"], mode="tag")
    )
    part.parent = SimpleNamespace(task_name="fdr")

    params = part._get_params_run(SimpleNamespace())

    assert (params["mode"], params["filter_name"]) == ("tag", "fdr")


def test_the_fdr_control_needs_varlociraptor_calls():
    config = VariantFiltration.model_construct(tool="varlociraptor")
    calls = DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"}))

    with pytest.raises(ValueError, match="produces variants tagged indel, snv, somatic"):
        VariantFiltrationWorkflow.task_produces(config, {"variants": calls})

    produced = VariantFiltrationWorkflow.task_produces(
        config, {"variants": calls.with_tags("varlociraptor")}
    )
    assert produced == (calls.with_tags("varlociraptor", "filtered"),)
    assert VariantFiltrationWorkflow.filter_ids(config) == frozenset()
