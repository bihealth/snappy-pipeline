# -*- coding: utf-8 -*-
"""Tests for step parts of ``somatic_neoepitope_prediction``."""

from types import SimpleNamespace

from snakemake.pathvars import Pathvars

from snappy_pipeline.workflows.somatic_neoepitope_prediction import NetChopStepPart
from snappy_pipeline.workflows.somatic_neoepitope_prediction.model import PvacTools


def test_container_is_an_image_uri_or_a_local_image_file():
    # The context that load_project() gives: a local image is a path of the configuration
    context = {"project_dir": "/projects/p1", "pathvars": Pathvars.with_defaults(), "paths": []}

    default = PvacTools.model_validate({}, context=context)
    local = PvacTools.model_validate({"container": "images/pvactools.sif"}, context=context)

    assert default.container == (
        "docker://griffithlab/pvactools"
        "@sha256:0297856ce7a9dec7339f3405f92fcfb37be48a6a2b06ec809270b93c7a81f25c"
    )
    assert local.container == "/projects/p1/images/pvactools.sif"
    assert context["paths"] == ["/projects/p1/images/pvactools.sif"]


def test_netchop_reads_the_options_of_its_prediction():
    # Several predictions run in one task; the rule's wildcard names the one to post-process
    part = object.__new__(NetChopStepPart)
    net_chop = SimpleNamespace(method="20s", threshold=0.7)
    part.parent = SimpleNamespace(
        sub_steps={"pvacfuse": SimpleNamespace(cfg=SimpleNamespace(net_chop=net_chop))}
    )

    params = part._get_params_run(SimpleNamespace(prediction="pvacfuse", tumor_library="T1"))

    assert params == {"tool": "pvacfuse", "method": "20s", "threshold": 0.7}
