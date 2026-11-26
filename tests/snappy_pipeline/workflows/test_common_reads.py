# -*- coding: utf-8 -*-
"""Tests for ``common/reads.py``: reads are inputs, unless another task writes them."""

from types import SimpleNamespace

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
        reads_input=lambda library: [f"tasks/{source}/work/{library}/out"],
        reads_from_task=lambda source: source_step != "external_data",
    )


def test_files_of_an_external_data_task_are_inputs():
    step = _step("raw", "external_data")
    assert reads_input_files(step, "lib") == {
        "reads_left": ["/raw/a_R1.fq.gz", "/raw/b_R1.fq.gz"],
        "reads_right": ["/raw/a_R2.fq.gz", "/raw/b_R2.fq.gz"],
    }
    assert reads_params(step, "lib") == {}


def test_single_end_has_no_right_reads():
    inputs = reads_input_files(_step("raw", "external_data", mates=("left",)), "lib")
    assert list(inputs) == ["reads_left"]


def test_reads_written_by_another_task_are_params():
    step = _step("trimming")
    assert reads_input_files(step, "lib") == {"reads": ["tasks/trimming/work/lib/out"]}
    assert reads_params(step, "lib")["input"]["reads_left"] == [
        "/raw/a_R1.fq.gz",
        "/raw/b_R1.fq.gz",
    ]
