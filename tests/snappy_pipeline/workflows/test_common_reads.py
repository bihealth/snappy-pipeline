# -*- coding: utf-8 -*-
"""Tests for ``common/reads.py``: reads are inputs, unless another task writes them."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params


def _step(source, source_step="adapter_trimming", mates=("left", "right")):
    groups = [
        SimpleNamespace(
            left="/raw/a_R1.fq.gz", right="/raw/a_R2.fq.gz" if "right" in mates else None
        ),
        SimpleNamespace(
            left="/raw/b_R1.fq.gz", right="/raw/b_R2.fq.gz" if "right" in mates else None
        ),
    ]
    return SimpleNamespace(
        depends_on=SimpleNamespace(reads=source),
        project=SimpleNamespace(task=lambda name: SimpleNamespace(step=source_step)),
        read_groups=lambda library: groups,
        reads_input=lambda library: [f"tasks/{source}/output/{library}/out/.done"],
    )


@pytest.mark.parametrize("source, step", [("data_sets", ""), ("raw", "external_data")])
def test_files_are_inputs(source, step):
    step = _step(source, step)
    assert reads_input_files(step, "lib") == {
        "reads_left": ["/raw/a_R1.fq.gz", "/raw/b_R1.fq.gz"],
        "reads_right": ["/raw/a_R2.fq.gz", "/raw/b_R2.fq.gz"],
    }
    assert reads_params(step, "lib") == {}


def test_single_end_has_no_right_reads():
    inputs = reads_input_files(_step("data_sets", mates=("left",)), "lib")
    assert list(inputs) == ["reads_left"]


def test_reads_written_by_another_task_are_params():
    step = _step("trimming")
    assert reads_input_files(step, "lib") == {"reads": ["tasks/trimming/output/lib/out/.done"]}
    assert reads_params(step, "lib")["input"]["reads_left"] == [
        "/raw/a_R1.fq.gz",
        "/raw/b_R1.fq.gz",
    ]
