# -*- coding: utf-8 -*-
"""Tests for the input functions of ``gene_expression_quantification`` step parts."""

from types import SimpleNamespace

from snappy_pipeline.workflows.gene_expression_quantification import (
    QCStepPartDupradar,
    QCStepPartRnaseqc,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

WILDCARDS = SimpleNamespace(library_name="L1")

ALIGNMENTS = {"bam": "tasks/mapping/output/L1.bam", "bai": "tasks/mapping/output/L1.bam.bai"}


def _part(cls, **tool_config):
    # The input functions only read upstream paths and config, so skip the full step setup.
    part = object.__new__(cls)
    part.parent = SimpleNamespace(
        get_upstream_paths=lambda field, library_name: ExpectedAlignments(**ALIGNMENTS)
    )
    part.config = SimpleNamespace(tool=cls.name, **{cls.name: SimpleNamespace(**tool_config)})
    part.w_config = SimpleNamespace(
        static_data_config=SimpleNamespace(reference=SimpleNamespace(path="/refs/genome.fa"))
    )
    return part


def test_dupradar_inputs_are_named():
    part = _part(QCStepPartDupradar, dupradar_path_annotation_gtf="genes.gtf")

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | {
        "dupradar_path_annotation_gtf": "genes.gtf"
    }


def test_rnaseqc_inputs_are_named():
    part = _part(QCStepPartRnaseqc, rnaseqc_path_annotation_gtf="genes.gtf")

    assert part._get_input_files_run(WILDCARDS) == ALIGNMENTS | {
        "reference": "/refs/genome.fa",
        "rnaseqc_path_annotation_gtf": "genes.gtf",
    }
