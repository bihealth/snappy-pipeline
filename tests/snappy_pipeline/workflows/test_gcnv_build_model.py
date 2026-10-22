# -*- coding: utf-8 -*-
"""Tests for the gCNV cohort-mode step parts in ``common/gcnv/gcnv_build_model.py``."""

from types import MethodType, SimpleNamespace

import pandas as pd
import pytest
from snakemake.io import is_flagged
from snakemake.iocontainers import Wildcards

from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import ExpectedReference
from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import (
    BuildGcnvModelStepPart,
    CallCnvsMixin,
    ContigPloidyMixin,
)
from snappy_pipeline.workflows.common.gcnv.gcnv_common import PreprocessIntervalsCommonMixin
from snappy_pipeline.workflows.gcnv_model import GcnvModelStepPart, GcnvModelWorkflow


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


def test_model_parts_are_directories_that_the_callers_read():
    """The ploidy model and the model of each shard are directories, not temp: other tasks
    read them through the contract."""
    step = object.__new__(BaseStep)
    step.config = SimpleNamespace(delete_intermediates=True)
    part = object.__new__(BuildGcnvModelStepPart)
    step._get_sub_step = lambda name: part

    outputs = step.get_output_files("gcnv", "contig_ploidy")
    outputs |= step.get_output_files("gcnv", "call_cnvs")

    assert outputs == {
        "ploidy": "work/{library_kit}/out/{library_kit}.contig_ploidy",
        "calls": "work/{library_kit}/out/{library_kit}.{shard}.call_cnvs",
    }
    assert all(is_flagged(path, "directory") for path in outputs.values())
    assert not any(is_flagged(path, "temp") for path in outputs.values())
    paths = GcnvModelWorkflow.get_output_paths(None, library_kit="kit", shard="0001_of_0001")
    assert paths["ploidy"] == "work/kit/out/kit.contig_ploidy"
    assert paths["calls"] == "work/kit/out/kit.0001_of_0001.call_cnvs"


def test_call_cnvs_reads_the_ploidy_directory():
    part = CallCnvsMixin()
    part.parent = SimpleNamespace(
        primary_dna_libraries=lambda include_background: pd.DataFrame({"library_name": ["P001"]})
    )
    part.ngs_library_to_kit = {"P001": "kit"}
    inputs = part._get_input_files_call_cnvs(SimpleNamespace(library_kit="kit"))
    assert inputs["ploidy"] == "work/kit/out/kit.contig_ploidy"


class _Checkpoint:
    def __init__(self, scatter_dir):
        self.scatter_dir = scatter_dir

    def get(self, **wildcards):
        return SimpleNamespace(output=[str(self.scatter_dir)])


def test_post_germline_calls_reads_the_model_directories(tmp_path):
    for shard in ("0001_of_0002", "0002_of_0002"):
        (tmp_path / f"temp_{shard}").mkdir()
        (tmp_path / f"temp_{shard}" / "scattered.interval_list").touch()
    part = object.__new__(GcnvModelStepPart)
    part.ngs_library_to_kit = {"P001": "kit"}
    checkpoints = SimpleNamespace(gcnv_model_scatter_intervals=_Checkpoint(tmp_path))

    inputs = part._get_input_files_post_germline_calls(
        Wildcards(fromdict={"library_name": "P001"}), checkpoints
    )

    assert sorted(inputs["calls"]) == [
        "work/kit/out/kit.0001_of_0002.call_cnvs",
        "work/kit/out/kit.0002_of_0002.call_cnvs",
    ]
    assert inputs["ploidy"] == "work/kit/out/kit.contig_ploidy"


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
