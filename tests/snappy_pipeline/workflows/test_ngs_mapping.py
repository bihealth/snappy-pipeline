# -*- coding: utf-8 -*-
"""Tests for the config model of the ``ngs_mapping`` step."""

from types import SimpleNamespace

import pandas as pd
import pydantic
import pytest

from snappy_pipeline.models import SnappyStepModel
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import NgsMapping

STAR_DEPENDS_ON = {
    "reads": "reads",
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


@pytest.mark.parametrize(
    "tool, mapped",
    [("bwa", ["dna"]), ("bwa_mem3", ["dna"]), ("minimap2", ["long"]), ("star", ["rna"])],
)
def test_each_mapper_takes_the_libraries_it_maps(tool, mapped):
    libraries = pd.DataFrame(
        [
            {"library_name": "dna", "extraction_type": "dna", "seq_platform": "Illumina"},
            {"library_name": "long", "extraction_type": "dna", "seq_platform": "PacBio"},
            {"library_name": "rna", "extraction_type": "rna", "seq_platform": "Illumina"},
        ]
    )
    requirement = NgsMappingWorkflow.task_library_requirement(SimpleNamespace(tool=tool))

    assert libraries.query(requirement)["library_name"].tolist() == mapped


def test_only_steps_with_a_bam_for_other_tasks_delete_it():
    class NoBam(SnappyStepModel):
        pass

    with pytest.raises(pydantic.ValidationError, match="writes no BAM file to delete"):
        NoBam(delete_intermediates=["bam"])
    assert NoBam(delete_intermediates=["intermediates"]).delete_intermediates == ["intermediates"]
    depends_on = {"reads": "reads", "reference": "genome", "index": "bwa"}
    config = NgsMapping(tool="bwa", bwa={}, depends_on=depends_on, delete_intermediates=["bam"])
    assert config.delete_intermediates == ["bam"]
