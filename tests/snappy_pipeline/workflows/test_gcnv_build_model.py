# -*- coding: utf-8 -*-
"""Tests for the gCNV cohort-mode step parts in ``common/gcnv/gcnv_build_model.py``."""

from types import MethodType, SimpleNamespace

import pandas as pd
import pytest

from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import ExpectedReference
from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import ContigPloidyMixin
from snappy_pipeline.workflows.common.gcnv.gcnv_common import PreprocessIntervalsCommonMixin


def _part(par_regions):
    part = ContigPloidyMixin()
    part.config = SimpleNamespace(depends_on=SimpleNamespace(par_regions=par_regions))
    part.parent = SimpleNamespace(
        get_upstream_paths=lambda field: SimpleNamespace(bed=f"/refs/{field}.bed"),
        primary_dna_libraries=lambda include_background: pd.DataFrame({"library_name": ["P001"]}),
        pedigrees=lambda include_background: {"P001": ["P001"]},
    )
    part.ngs_library_to_kit = {"P001": "kit"}
    return part


@pytest.mark.parametrize("path, expected", [("par", "/refs/par_regions.bed"), ("", None)])
def test_contig_ploidy_reads_par_intervals_as_input(path, expected):
    inputs = _part(path)._get_input_files_contig_ploidy(SimpleNamespace(library_kit="kit"))
    assert inputs.get("par_intervals") == expected
    assert inputs["tsv"] == ["work/P001/out/P001.coverage.tsv"]


#: The genome and its sidecar files as rule inputs
REFERENCE_INPUTS = {
    "reference": "/refs/genome.fa",
    "reference_fai": "/refs/genome.fa.fai",
    "reference_dict": "/refs/genome.dict",
}


def _preprocess_part(intervals):
    part = PreprocessIntervalsCommonMixin()
    paths = {
        "reference": ExpectedReference(fasta="/refs/genome.fa"),
        "kits": SimpleNamespace(targets="/refs/kit.bed"),
    }
    part.parent = SimpleNamespace(get_upstream_paths=lambda field, **kwargs: paths[field])
    part.parent.reference_inputs = MethodType(BaseStep.reference_inputs, part.parent)
    part.config = SimpleNamespace(intervals=intervals)
    return part


def test_preprocess_intervals_reads_the_target_bed_of_the_kit():
    inputs = _preprocess_part("targeted")._get_input_files_preprocess_intervals(
        SimpleNamespace(library_kit="kit")
    )
    assert inputs == REFERENCE_INPUTS | {"target_bed": "/refs/kit.bed"}


def test_preprocess_intervals_without_targets_reads_only_the_reference():
    inputs = _preprocess_part("wgs")._get_input_files_preprocess_intervals(
        SimpleNamespace(library_kit="wgs")
    )
    assert inputs == REFERENCE_INPUTS
