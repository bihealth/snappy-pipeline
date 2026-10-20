# -*- coding: utf-8 -*-
"""Tests for ``BaseStep.read_groups``: data sets and ``external_data`` tasks find reads the same way."""

from types import SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.reads import ReadDiscovery
from snappy_pipeline.workflow_model import DataSet
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.external_data.model import ExternalData

PATTERN = {"left": r"(?P<readgroup>.+)_R1\.fq\.gz", "right": r"(?P<readgroup>.+)_R2\.fq\.gz"}
RAW = {"type": "raw", "tags": []}


def _touch(root, *names):
    for name in names:
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).touch()


def _step(source, data_sets=None, tasks=None):
    # read_groups only reads these attributes, so skip the full step setup
    step = object.__new__(BaseStep)
    step.depends_on = SimpleNamespace(reads=source)
    step.w_config = SimpleNamespace(data_sets=data_sets or {})
    step.data_set_infos = [SimpleNamespace(name=name) for name in data_sets or {}]
    step.project = SimpleNamespace(
        read_discovery=ReadDiscovery(),
        task=lambda name: SimpleNamespace(step="external_data"),
        task_configs=tasks or {},
    )
    return step


def test_data_set_reads_are_found_only_below_the_library_data_set(tmp_path):
    _touch(tmp_path / "batch1", "lib/a_R1.fq.gz", "lib/a_R2.fq.gz")
    _touch(tmp_path / "batch2", "lib/b_R1.fq.gz", "lib/b_R2.fq.gz")
    data_sets = {
        name: DataSet(search_paths=[str(tmp_path / name)], search_patterns=[PATTERN])
        for name in ("batch1", "batch2")
    }
    step = _step("data_sets", data_sets=data_sets)
    # The library is in the sheet of batch2, in the folder "lib"
    step.__dict__["_library_rows"] = {"L1": [{"data_set": "batch2", "folder_name": "lib"}]}
    groups = step.read_groups("L1")
    assert [group.left for group in groups] == [str(tmp_path / "batch2/lib/b_R1.fq.gz")]


def test_external_reads_may_mix_single_and_paired_end(tmp_path):
    _touch(tmp_path, "L1/a_R1.fq.gz")
    config = {"search_paths": [str(tmp_path)], "search_patterns": [PATTERN], "produces": RAW}
    with pytest.raises(ValueError, match="set mixed_se_pe"):
        _step("raw", tasks={"raw": ExternalData(**config)}).read_groups("L1")

    tasks = {"raw": ExternalData(**config, mixed_se_pe=True)}
    groups = _step("raw", tasks=tasks).read_groups("L1")
    assert [(group.left, group.right) for group in groups] == [
        (str(tmp_path / "L1/a_R1.fq.gz"), None)
    ]


def test_mixed_se_pe_is_only_for_reads():
    config = {"search_paths": ["/vcfs"], "search_patterns": [{"vcf": r".+\.vcf\.gz"}]}
    with pytest.raises(pydantic.ValidationError, match="mixed_se_pe is only for reads"):
        ExternalData(**config, produces={"type": "variants", "tags": []}, mixed_se_pe=True)
