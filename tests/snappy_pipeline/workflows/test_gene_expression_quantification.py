# -*- coding: utf-8 -*-
"""Tests for the input functions of ``gene_expression_quantification`` step parts."""

import typing
from types import MethodType, SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.workflows.gene_expression_quantification import (
    RNASEQC_FEATURES_SIGNATURE,
    QCStepPartDupradar,
    QCStepPartRnaseqc,
)
from snappy_pipeline.workflows.gene_expression_quantification.model import (
    ExpectedStrandedness,
    GeneExpressionQuantification,
)
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import ExpectedFeatures, ExpectedReference
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

WILDCARDS = SimpleNamespace(library_name="L1")

ALIGNMENTS = {"bam": "tasks/mapping/output/L1.bam", "bai": "tasks/mapping/output/L1.bam.bai"}

DECISION = {"decision": "tasks/strandedness/output/L1/out/L1.decision"}

UPSTREAM = {
    "alignments": ExpectedAlignments(**ALIGNMENTS),
    "strandedness": ExpectedStrandedness(**DECISION),
    "reference": ExpectedReference(fasta="/refs/genome.fa"),
    "features": ExpectedFeatures(gtf="/refs/genes.gtf"),
}


def _part(cls, **tool_config):
    # The input functions only read upstream paths and config, so skip the full step setup.
    part = object.__new__(cls)
    signatures = {}

    def get_upstream_paths(field, signature=None, **kwargs):
        signatures[field] = signature
        return UPSTREAM[field]

    part.parent = SimpleNamespace(get_upstream_paths=get_upstream_paths, signatures=signatures)
    part.parent.reference_inputs = MethodType(BaseStep.reference_inputs, part.parent)
    part.config = SimpleNamespace(tool=cls.name, **{cls.name: SimpleNamespace(**tool_config)})
    return part


def test_dupradar_inputs_are_named():
    part = _part(QCStepPartDupradar)

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | DECISION | {
        "features": "/refs/genes.gtf"
    }


def test_rnaseqc_inputs_are_named():
    part = _part(QCStepPartRnaseqc)

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | DECISION | {
        "reference": "/refs/genome.fa",
        "reference_fai": "/refs/genome.fa.fai",
        "reference_dict": "/refs/genome.dict",
        "features": "/refs/genes.gtf",
    }
    # The GTF must be one prepared for RNA-SeQC
    assert part.parent.signatures["features"] == RNASEQC_FEATURES_SIGNATURE


def _tool_section(tool):
    """Return a minimal config section for ``tool``: placeholders for its required fields."""
    annotation = GeneExpressionQuantification.model_fields[tool].annotation
    section_model = next(a for a in typing.get_args(annotation) if a is not type(None))
    return {
        name: "placeholder"
        for name, field in section_model.model_fields.items()
        if field.is_required()
    }


@pytest.mark.parametrize("tool", ["featurecounts", "dupradar", "duplication", "rnaseqc", "stats"])
def test_quantifiers_require_a_strandedness_task(tool):
    section = {tool: _tool_section(tool)}
    with pytest.raises(pydantic.ValidationError, match="needs depends_on.strandedness"):
        GeneExpressionQuantification(tool=tool, depends_on={"alignments": "mapping"}, **section)

    depends_on = {"alignments": "mapping", "strandedness": "strandedness"}
    depends_on |= {
        "featurecounts": {"features": "genes"},
        "dupradar": {"features": "genes"},
        "rnaseqc": {"reference": "genome", "features": "genes"},
    }.get(tool, {})
    config = GeneExpressionQuantification(tool=tool, depends_on=depends_on, **section)
    assert config.depends_on.strandedness == "strandedness"


def test_salmon_requires_reads():
    with pytest.raises(pydantic.ValidationError, match="tool=salmon needs depends_on.reads"):
        GeneExpressionQuantification(tool="salmon", salmon=_tool_section("salmon"))


def test_strandedness_requires_exons():
    depends_on = {"alignments": "mapping"}
    with pytest.raises(pydantic.ValidationError, match="tool=strandedness needs depends_on.exons"):
        GeneExpressionQuantification(tool="strandedness", depends_on=depends_on, strandedness={})
