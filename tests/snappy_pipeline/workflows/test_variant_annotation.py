"""Tests for the VEP summary output of ``variant_annotation``."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.models.annotation import Vep
from snappy_pipeline.workflows.variant_annotation import VepStepPart


def _outputs(**cfg):
    part = object.__new__(VepStepPart)
    part.config = SimpleNamespace(vep=Vep(**cfg))
    return part._get_output_files_run()


def test_summary_is_an_output_by_default():
    outputs = _outputs()
    assert outputs["summary"] == "work/{library_name}/report/{library_name}.summary.html"
    assert "summary" in VepStepPart.report_outputs["run"]


def test_no_stats_declares_no_summary():
    assert "summary" not in _outputs(output_options=["everything", "no_stats"])


def test_stats_flags_are_not_extra_arguments():
    with pytest.raises(ValueError):
        Vep(extra_arguments=["--no_stats"])
