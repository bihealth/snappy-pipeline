# -*- coding: utf-8 -*-
"""Tests for step parts of ``somatic_neoepitope_prediction``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_neoepitope_prediction import NetChopStepPart


def test_netchop_params_use_the_task_tool():
    # Rule wildcards hold no tool; the tool is the task's configured one.
    part = object.__new__(NetChopStepPart)
    net_chop = SimpleNamespace(method="cterm", threshold=0.5)
    part.config = SimpleNamespace(
        tool="pvacseq", get=lambda name: SimpleNamespace(net_chop=net_chop)
    )

    params = part._get_params_pvacseq(SimpleNamespace(tumor_dna="T1"))

    assert params == {"tool": "pvacseq", "method": "cterm", "threshold": 0.5}
