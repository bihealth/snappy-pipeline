# -*- coding: utf-8 -*-
"""Tests for ``BaseStep.read_groups`` with the files tables of an ``external_data`` task."""

from types import SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.base import ConfigError
from snappy_pipeline.reads import read_files_tables
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.external_data.model import ExternalData

RAW = {"type": "raw", "tags": []}


def _step(source, tasks):
    # read_groups only reads these attributes, so skip the full step setup
    step = object.__new__(BaseStep)
    step.depends_on = SimpleNamespace(reads=source)
    step.project = SimpleNamespace(
        task=lambda name: SimpleNamespace(step="external_data"),
        task_configs=tasks,
        files_table=lambda config: read_files_tables(
            config.library_files, True, config.mixed_se_pe
        ),
    )
    return step


def test_external_reads_may_mix_single_and_paired_end(tmp_path):
    table = tmp_path / "reads.tsv"
    table.write_text(
        "library\treadgroup\tleft\tright\nL1\ta\ta_R1.fq.gz\t\nL2\tb\tb_R1.fq.gz\tb_R2.fq.gz\n"
    )
    config = {"library_files": [str(table)], "produces": RAW}
    with pytest.raises(ValueError, match="set mixed_se_pe"):
        _step("raw", {"raw": ExternalData(**config)}).read_groups("L1")

    tasks = {"raw": ExternalData(**config, mixed_se_pe=True)}
    groups = _step("raw", tasks).read_groups("L1")
    assert [(group.left, group.right) for group in groups] == [(str(tmp_path / "a_R1.fq.gz"), None)]
    with pytest.raises(ConfigError, match="list no files of library 'L3'; run `snappy scan-files`"):
        _step("raw", tasks).read_groups("L3")


def test_mixed_se_pe_is_only_for_reads():
    config = {"library_files": ["vcfs.tsv"], "produces": {"type": "variants", "tags": []}}
    with pytest.raises(pydantic.ValidationError, match="mixed_se_pe is only for reads"):
        ExternalData(**config, mixed_se_pe=True)
