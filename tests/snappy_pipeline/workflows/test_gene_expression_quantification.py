# -*- coding: utf-8 -*-
"""Tests for the input functions of ``gene_expression_quantification`` step parts."""

import typing
from types import SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.workflows.gene_expression_quantification import (
    QCStepPartDupradar,
    QCStepPartRnaseqc,
)
from snappy_pipeline.workflows.gene_expression_quantification.model import (
    ExpectedStrandedness,
    GeneExpressionQuantification,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

WILDCARDS = SimpleNamespace(library_name="L1")

ALIGNMENTS = {"bam": "tasks/mapping/output/L1.bam", "bai": "tasks/mapping/output/L1.bam.bai"}

DECISION = {"decision": "tasks/strandedness/output/L1/out/L1.decision"}

UPSTREAM = {
    "alignments": ExpectedAlignments(**ALIGNMENTS),
    "strandedness": ExpectedStrandedness(**DECISION),
}


def _part(cls, **tool_config):
    # The input functions only read upstream paths and config, so skip the full step setup.
    part = object.__new__(cls)
    part.parent = SimpleNamespace(get_upstream_paths=lambda field, library_name: UPSTREAM[field])
    part.config = SimpleNamespace(tool=cls.name, **{cls.name: SimpleNamespace(**tool_config)})
    part.w_config = SimpleNamespace(
        static_data_config=SimpleNamespace(reference=SimpleNamespace(path="/refs/genome.fa"))
    )
    return part


def test_dupradar_inputs_are_named():
    part = _part(QCStepPartDupradar, dupradar_path_annotation_gtf="genes.gtf")

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | DECISION | {
        "dupradar_path_annotation_gtf": "genes.gtf"
    }


def test_rnaseqc_inputs_are_named():
    part = _part(QCStepPartRnaseqc, rnaseqc_path_annotation_gtf="genes.gtf")

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | DECISION | {
        "reference": "/refs/genome.fa",
        "rnaseqc_path_annotation_gtf": "genes.gtf",
    }


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
    with pytest.raises(pydantic.ValidationError, match="needs depends_on.strandedness"):
        GeneExpressionQuantification(tool=tool, **{tool: _tool_section(tool)})

    config = GeneExpressionQuantification(
        tool=tool, depends_on={"strandedness": "strandedness"}, **{tool: _tool_section(tool)}
    )
    assert config.depends_on.strandedness == "strandedness"
