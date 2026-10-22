# -*- coding: utf-8 -*-
"""Tests for the tool settings that ``ngs_mapping`` passes to its wrappers."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows import ngs_mapping
from snappy_pipeline.workflows.ngs_mapping import BwaMem2StepPart, BwaStepPart
from snappy_pipeline.workflows.ngs_mapping.model import Bwa, BwaMem2


@pytest.fixture(autouse=True)
def no_reads(monkeypatch):
    """The read params need the whole step; the settings do not."""
    monkeypatch.setattr(ngs_mapping, "reads_params", lambda parent, library_name: {})


def _params(part_class, **config):
    part = object.__new__(part_class)
    part.config = SimpleNamespace(**config)
    part.parent = None
    return part._get_params_run(SimpleNamespace(library_name="lib"))


def test_bwa_passes_platform_and_aln_threshold():
    params = _params(BwaStepPart, bwa=Bwa(platform="ONT", aln_max_read_length=50))

    assert params["platform"] == "ONT"
    assert params["aln_max_read_length"] == 50


def test_bwa_defaults_keep_the_former_literals():
    params = _params(BwaStepPart, bwa=Bwa())

    assert params["platform"] == "ILLUMINA"
    assert params["aln_max_read_length"] == 75


def test_bwa_mem2_passes_platform():
    params = _params(BwaMem2StepPart, bwa_mem2=BwaMem2(platform="ONT"))

    assert params["platform"] == "ONT"
