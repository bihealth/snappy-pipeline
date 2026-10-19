# -*- coding: utf-8 -*-
"""Tests for the gCNV case mode in ``common/gcnv/gcnv_run.py`` and the kit mapping."""

from types import SimpleNamespace

import pandas as pd
import pytest

from snappy_pipeline.workflows.common.gcnv.gcnv_common import build_ngs_library_to_kit
from snappy_pipeline.workflows.common.gcnv.gcnv_run import (
    CallCnvsMixin,
    PostGermlineCallsMixin,
    ShardsMixin,
)
from snappy_pipeline.workflows.gcnv_model.model import GcnvModel


def _upstream(field, library_kit, shard="{shard}"):
    assert field == "gcnv_model"
    prefix = f"tasks/model/work/{library_kit}/out/{library_kit}"
    return SimpleNamespace(calls=f"{prefix}.{shard}.call_cnvs/.done")


class _Checkpoint:
    def __init__(self, shards_file):
        self.shards_file = shards_file

    def get(self, **wildcards):
        assert wildcards == {"library_kit": "wgs"}
        return SimpleNamespace(output=[str(self.shards_file)])


def _part(cls):
    part = cls()
    part.parent = SimpleNamespace(name="sv_calling_wgs", get_upstream_paths=_upstream)
    part.ngs_library_to_kit = {"P001": "wgs"}
    return part


def test_shards_are_listed_from_the_model_directory(tmp_path):
    for name in ("0002_of_0002", "0001_of_0002"):
        (tmp_path / "scatter" / f"temp_{name}").mkdir(parents=True)
    (tmp_path / "scatter" / "file.txt").touch()
    out = tmp_path / "shards.txt"
    ShardsMixin().write_shards(
        SimpleNamespace(shards=str(tmp_path / "scatter")), SimpleNamespace(shards=str(out))
    )
    assert out.read_text() == "0001_of_0002\n0002_of_0002\n"


def test_post_germline_calls_reads_the_shards_from_the_local_checkpoint(tmp_path):
    shards_file = tmp_path / "shards.txt"
    shards_file.write_text("0001_of_0002\n0002_of_0002\n")
    checkpoints = SimpleNamespace(sv_calling_wgs_gcnv_shards=_Checkpoint(shards_file))
    inputs = _part(PostGermlineCallsMixin)._get_input_files_post_germline_calls(
        SimpleNamespace(library_name="P001"), checkpoints
    )
    calls = [f"work/wgs/out/wgs.{s}.call_cnvs/.done" for s in ("0001_of_0002", "0002_of_0002")]
    assert inputs["calls"] == calls
    assert inputs["ploidy"] == "work/wgs/out/wgs.contig_ploidy/.done"
    assert inputs["model"] == [f"tasks/model/{c}" for c in calls]


def test_call_cnvs_reads_the_model_of_its_shard():
    part = _part(CallCnvsMixin)
    part.parent.primary_dna_libraries = lambda include_background: pd.DataFrame(
        {"library_name": ["P001"]}
    )
    inputs = part._get_input_files_call_cnvs(
        SimpleNamespace(library_kit="wgs", shard="0001_of_0002")
    )
    assert inputs["model"] == "tasks/model/work/wgs/out/wgs.0001_of_0002.call_cnvs/.done"
    assert inputs["tsv"] == ["work/P001/out/P001.coverage.tsv"]


def _model_config(intervals, kits=""):
    return GcnvModel.model_validate(
        {
            "depends_on": {
                "alignments": "mapping",
                "reference": "genome",
                "kits": kits,
                "mappable_regions": "mappable",
            },
            "intervals": intervals,
        }
    )


LIBRARIES = pd.DataFrame({"library_name": ["A", "B", "C"]})
LIBRARY_KITS = {"A": "Kit_v1", "B": "Kit_v2"}


def test_kit_mapping_of_targeted_models_is_the_library_kit():
    assert build_ngs_library_to_kit(LIBRARIES, "targeted", LIBRARY_KITS) == {
        "A": "Kit_v1",
        "B": "Kit_v2",
    }


def test_kit_of_wgs_models_is_wgs():
    assert build_ngs_library_to_kit(LIBRARIES, "wgs", LIBRARY_KITS) == {
        "A": "wgs",
        "B": "wgs",
        "C": "wgs",
    }


@pytest.mark.parametrize("intervals, kits", [("targeted", ""), ("wgs", "kits")])
def test_model_config_needs_kits_only_for_targeted_intervals(intervals, kits):
    with pytest.raises(ValueError, match="depends_on.kits"):
        _model_config(intervals, kits)
