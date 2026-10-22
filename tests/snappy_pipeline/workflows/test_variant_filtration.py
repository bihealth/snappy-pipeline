# -*- coding: utf-8 -*-
"""Tests for the params of ``variant_filtration``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.variant_filtration import (
    EbfilterStepPart,
    VariantFiltrationStepPart,
    VembraneStepPart,
)
from snappy_pipeline.workflows.variant_filtration.model import Ebfilter, Vembrane


def _part(cls, tool, section, annotated=True):
    # The params only read the tool's config section, so skip the full step setup.
    part = object.__new__(cls)
    part.config = SimpleNamespace(tool=tool, **{tool: section})
    part.parent = SimpleNamespace(task_name="filter_task", variants_are_annotated=annotated)
    return part


@pytest.mark.parametrize("annotated", [True, False])
def test_ebfilter_knows_whether_the_variants_are_annotated(annotated):
    part = _part(EbfilterStepPart, "ebfilter", Ebfilter(), annotated)

    params = part._get_params_run(SimpleNamespace())

    assert params["has_annotation"] is annotated
    assert params["filter_name"] == "filter_task"


def test_ebfilter_passes_what_skips_it():
    part = _part(EbfilterStepPart, "ebfilter", Ebfilter(skip_annotation="intron_variant"))

    params = part._get_params_run(SimpleNamespace())

    assert params["skip_filters"] == ["germline", "weak_evidence", "OffExome"]
    assert params["skip_annotation"] == "intron_variant"


@pytest.mark.parametrize("mode", ["tag", "filter"])
def test_vembrane_passes_aux_context_and_ontology_in_both_modes(mode, monkeypatch):
    expressions = (
        {"expressions": {"low": "QUAL < 30"}} if mode == "tag" else {"expression": "QUAL >= 30"}
    )
    section = Vembrane(
        mode=mode,
        aux={"genes": "genes.txt"},
        context=["MIN_DP = 10"],
        context_files=["helpers.py"],
        ontology="so.obo",
        **expressions,
    )
    part = _part(VembraneStepPart, "vembrane", section)
    monkeypatch.setattr(
        VariantFiltrationStepPart, "_get_input_files_run", lambda self, wildcards: {}
    )

    params = part._get_params_run(SimpleNamespace())
    inputs = part._get_input_files_run(SimpleNamespace())

    assert (params["aux_names"], params["context"]) == (["genes"], ["MIN_DP = 10"])
    assert set(inputs) == {"aux", "context_files", "ontology"}
