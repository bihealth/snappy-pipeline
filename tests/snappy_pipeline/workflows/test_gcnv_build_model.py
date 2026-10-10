# -*- coding: utf-8 -*-
"""Tests for the gCNV cohort-mode step parts in ``common/gcnv/gcnv_build_model.py``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import ContigPloidyMixin


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
