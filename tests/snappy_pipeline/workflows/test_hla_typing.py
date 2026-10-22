# -*- coding: utf-8 -*-
"""Tests for the arcasHLA params of ``hla_typing``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.hla_typing import ArcasHlaStepPart
from snappy_pipeline.workflows.hla_typing.model import ArcasHla


def _arcashla_params(**cfg):
    # The params only read the arcasHLA config section, so skip the full step setup.
    part = object.__new__(ArcasHlaStepPart)
    part.config = SimpleNamespace(arcashla=ArcasHla(**cfg))
    return part._get_params_run(SimpleNamespace(library_name="lib"))


def test_arcashla_passes_its_genotype_options():
    assert _arcashla_params(population="caucasian") == {
        "population": "caucasian",
        "min_count": 75,
        "tolerance": 1e-7,
        "max_iterations": 1000,
        "drop_threshold": 0.1,
        "zygosity_threshold": 0.15,
    }


def test_arcashla_passes_unset_options_only_when_set():
    params = _arcashla_params(drop_iterations=10, avg=180, std=25)

    assert (params["drop_iterations"], params["avg"], params["std"]) == (10, 180, 25)


def test_optitype_passes_its_settings(monkeypatch):
    from snappy_pipeline.workflows import hla_typing
    from snappy_pipeline.workflows.hla_typing import OptiTypeStepPart
    from snappy_pipeline.workflows.hla_typing.model import Optitype

    monkeypatch.setattr(hla_typing, "reads_params", lambda parent, library_name: {})
    part = object.__new__(OptiTypeStepPart)
    part.parent = SimpleNamespace(extraction_types={"lib": "dna"})
    part.config = SimpleNamespace(optitype=Optitype(solver="cbc", unpaired_weight=0.5))

    params = part._get_params_run(SimpleNamespace(library_name="lib"))

    assert (params["solver"], params["unpaired_weight"]) == ("cbc", 0.5)
    assert Optitype().solver == "glpk" and Optitype().unpaired_weight == 0
