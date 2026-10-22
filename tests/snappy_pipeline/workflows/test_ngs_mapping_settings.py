# -*- coding: utf-8 -*-
"""Tests for the tool settings that ``ngs_mapping`` passes to its wrappers."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows import ngs_mapping
from snappy_pipeline.workflows.ngs_mapping import (
    BwaMem2StepPart,
    BwaStepPart,
    Minimap2StepPart,
    StarStepPart,
)
from snappy_pipeline.workflows.ngs_mapping.model import Bwa, BwaMem2, Minimap2, Star


@pytest.fixture(autouse=True)
def no_reads(monkeypatch):
    """The read params need the whole step; the settings do not."""
    monkeypatch.setattr(ngs_mapping, "reads_params", lambda parent, library_name: {})


def _params(part_class, **config):
    part = object.__new__(part_class)
    part.config = SimpleNamespace(**config)
    part.parent = SimpleNamespace(library_infos={"lib": {"library_kit": "PacBio HiFi"}})
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


def test_star_passes_platform():
    assert _params(StarStepPart, star=Star())["platform"] == "ILLUMINA"
    assert _params(StarStepPart, star=Star(platform="ONT"))["platform"] == "ONT"


def test_minimap2_defaults_keep_the_former_literals():
    params = _params(Minimap2StepPart, minimap2=Minimap2())

    assert params["platform"] == "PACBIO"
    assert params["soft_clip_supplementary"] is True
    assert params["output_md_tag"] is True
    assert params["memory_bam_sort"] == "4G"


def test_minimap2_passes_its_settings():
    params = _params(
        Minimap2StepPart,
        minimap2=Minimap2(
            platform="ONT",
            soft_clip_supplementary=False,
            output_md_tag=False,
            memory_bam_sort="2G",
        ),
    )

    assert params["platform"] == "ONT"
    assert params["soft_clip_supplementary"] is False
    assert params["output_md_tag"] is False
    assert params["memory_bam_sort"] == "2G"
