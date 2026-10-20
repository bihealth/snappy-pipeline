# -*- coding: utf-8 -*-
"""Tests for step parts of ``somatic_neoepitope_prediction``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_neoepitope_prediction import NetChopStepPart


def test_netchop_reads_the_options_of_its_prediction():
    # Several predictions run in one task; the rule's wildcard names the one to post-process
    part = object.__new__(NetChopStepPart)
    net_chop = SimpleNamespace(method="20s", threshold=0.7)
    part.parent = SimpleNamespace(
        sub_steps={"pvacfuse": SimpleNamespace(cfg=SimpleNamespace(net_chop=net_chop))}
    )

    params = part._get_params_run(SimpleNamespace(prediction="pvacfuse", tumor_library="T1"))

    assert params == {"tool": "pvacfuse", "method": "20s", "threshold": 0.7}
