# -*- coding: utf-8 -*-
"""Tests for the gCNV cohort-mode step parts in ``common/gcnv/gcnv_build_model.py``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import ContigPloidyMixin
from snappy_pipeline.workflows.common.gcnv.gcnv_common import PreprocessIntervalsCommonMixin


def _part(path_par_intervals):
    part = ContigPloidyMixin()
    part.config = SimpleNamespace(gcnv=SimpleNamespace(path_par_intervals=path_par_intervals))
    part.index_ngs_library_to_donor = {"P001": None}
    part.index_ngs_library_to_pedigree = {"P001": None}
    part.ngs_library_to_kit = {"P001": "kit"}
    return part


@pytest.mark.parametrize("path, expected", [("/refs/par.bed", "/refs/par.bed"), ("", None)])
def test_contig_ploidy_reads_par_intervals_as_input(path, expected):
    inputs = _part(path)._get_input_files_contig_ploidy(SimpleNamespace(library_kit="kit"))
    assert inputs.get("par_intervals") == expected
    assert inputs["tsv"] == ["work/P001/out/P001.coverage.tsv"]


def _preprocess_part(mapping):
    part = PreprocessIntervalsCommonMixin()
    reference = SimpleNamespace(fasta="/refs/genome.fa")
    part.parent = SimpleNamespace(get_upstream_paths=lambda field: reference)
    part.config = SimpleNamespace(gcnv=SimpleNamespace(path_target_interval_list_mapping=mapping))
    return part


def test_preprocess_intervals_reads_the_target_bed_of_the_kit():
    mapping = [SimpleNamespace(name="kit", path="/refs/kit.bed")]
    inputs = _preprocess_part(mapping)._get_input_files_preprocess_intervals(
        SimpleNamespace(library_kit="kit")
    )
    assert inputs == {"reference": "/refs/genome.fa", "target_bed": "/refs/kit.bed"}


def test_preprocess_intervals_without_targets_reads_only_the_reference():
    inputs = _preprocess_part(None)._get_input_files_preprocess_intervals(
        SimpleNamespace(library_kit="kit")
    )
    assert inputs == {"reference": "/refs/genome.fa"}
