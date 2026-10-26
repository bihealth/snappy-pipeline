# -*- coding: utf-8 -*-
"""Tests for the tool settings that ``ngs_data_qc`` passes to its wrappers."""

from types import SimpleNamespace

from snappy_pipeline.workflows.ngs_data_qc import AlfredStepPart
from snappy_pipeline.workflows.ngs_data_qc.model import Alfred


def test_alfred_passes_ignore_read_groups():
    for value in (True, False):
        part = object.__new__(AlfredStepPart)
        part.config = SimpleNamespace(alfred=Alfred(ignore_read_groups=value))

        assert part._get_params_run(SimpleNamespace(library_name="lib")) == {
            "ignore_read_groups": value
        }
