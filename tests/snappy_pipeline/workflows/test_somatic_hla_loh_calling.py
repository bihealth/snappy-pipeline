# -*- coding: utf-8 -*-
"""Tests for the LOHHLA params and outputs of ``somatic_hla_loh_calling``."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.workflows.somatic_hla_loh_calling import (
    LohhlaStepPart,
    SomaticHlaLohCallingWorkflow,
)


def test_lohhla_passes_the_minimum_coverage():
    # The params only read the step config, so skip the full step setup.
    part = object.__new__(LohhlaStepPart)
    part.config = SimpleNamespace(min_coverage_filter=20)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WGS1"))

    assert params == {"tumor_library": "T1-DNA1-WGS1", "min_coverage_filter": 20}


def test_lohhla_writes_its_output_directory():
    outputs = object.__new__(LohhlaStepPart)._get_output_files_run()
    assert outputs == {"out": "work/{tumor_library}/out"}
    assert is_flagged(outputs["out"], "directory")
    assert SomaticHlaLohCallingWorkflow.get_output_paths(None, library_name="T1") == {
        "out": "work/T1/out"
    }
