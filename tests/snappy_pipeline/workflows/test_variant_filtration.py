# -*- coding: utf-8 -*-
"""Tests for the params of ``variant_filtration``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.variant_filtration import EbfilterStepPart
from snappy_pipeline.workflows.variant_filtration.model import Ebfilter


def _part(cls, tool, section, annotated=True):
    # The params only read the tool's config section, so skip the full step setup.
    part = object.__new__(cls)
    part.config = SimpleNamespace(tool=tool, **{tool: section})
    part.parent = SimpleNamespace(task_name="filter_task", variants_are_annotated=annotated)
    return part


@pytest.mark.parametrize("annotated", [True, False])
def test_ebfilter_knows_whether_the_variants_are_annotated(annotated):
    part = _part(EbfilterStepPart, "ebfilter", Ebfilter(), annotated)

    params = part._get_params_run(SimpleNamespace())

    assert params["has_annotation"] is annotated
    assert params["filter_name"] == "filter_task"
