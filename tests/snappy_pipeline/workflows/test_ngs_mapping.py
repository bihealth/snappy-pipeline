# -*- coding: utf-8 -*-
"""Tests for the config model of the ``ngs_mapping`` step."""

import pydantic
import pytest

from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.workflows.ngs_mapping.model import NgsMapping

STAR_DEPENDS_ON = {
    "reads": "data_sets",
    "reference": "genome",
    "index": "star",
    "features": "genes",
}


def test_star_requires_exons():
    with pytest.raises(pydantic.ValidationError, match="tool=star needs depends_on.exons"):
        NgsMapping(tool="star", star={}, depends_on=STAR_DEPENDS_ON)


def test_star_infers_the_strand_without_a_strandedness_section():
    # The strandedness rule runs for every RNA library and its wrapper reads both values.
    config = NgsMapping(tool="star", star={}, depends_on=STAR_DEPENDS_ON | {"exons": "exons"})
    assert config.strandedness.model_dump() == {
        "strand": -1,
        "threshold": 0.85,
        "extra_arguments": [],
    }


MBCS_DEPENDS_ON = {"reads": "data_sets", "reference": "genome", "index": "bwa"}


def test_mbcs_recalibration_requires_common_variants():
    with pytest.raises(pydantic.ValidationError, match="needs depends_on.common_variants"):
        NgsMapping(
            tool="mbcs",
            mbcs={"mapping_tool": "bwa"},
            bwa={},
            depends_on=MBCS_DEPENDS_ON,
        )


def test_mbcs_without_recalibration_needs_no_common_variants():
    NgsMapping(
        tool="mbcs",
        mbcs={"mapping_tool": "bwa", "recalibrate": False},
        bwa={},
        depends_on=MBCS_DEPENDS_ON,
    )


def test_only_steps_with_a_bam_for_other_tasks_delete_it():
    class NoBam(SnappyStepModel):
        pass

    with pytest.raises(pydantic.ValidationError, match="writes no BAM file to delete"):
        NoBam(delete_intermediates=["bam"])
    assert NoBam(delete_intermediates=["intermediates"]).delete_intermediates == ["intermediates"]
    depends_on = {"reads": "data_sets", "reference": "genome", "index": "bwa"}
    config = NgsMapping(tool="bwa", bwa={}, depends_on=depends_on, delete_intermediates=["bam"])
    assert config.delete_intermediates == ["bam"]
