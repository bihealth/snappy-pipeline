# -*- coding: utf-8 -*-
"""Tests for the outputs of ``somatic_gene_fusion_calling``."""

from types import SimpleNamespace

import pytest
from snakemake.io import is_flagged

from snappy_pipeline.workflows.somatic_gene_fusion_calling import (
    ArribaStepPart,
    SomaticGeneFusionCallingWorkflow,
    StarFusionStepPart,
)
from snappy_pipeline.workflows.somatic_gene_fusion_calling.model import Tool


def test_star_fusion_writes_its_output_directory():
    outputs = object.__new__(StarFusionStepPart)._get_output_files_run()
    assert outputs == {"out": "work/{library_name}/out"}
    assert is_flagged(outputs["out"], "directory")


def test_arriba_writes_its_tables_without_a_marker():
    outputs = object.__new__(ArribaStepPart)._get_output_files_run()
    assert sorted(outputs) == ["discarded", "discarded_md5", "fusions", "fusions_md5"]


@pytest.mark.parametrize(
    "tool, expected",
    [
        (Tool.arriba, {"fusions": "work/L1/out/L1.fusions.tsv"}),
        (Tool.star_fusion, {"out": "work/L1/out"}),
    ],
)
def test_contract_names_the_outputs_of_the_tool(tool, expected):
    config = SimpleNamespace(tool=tool)
    paths = SomaticGeneFusionCallingWorkflow.get_output_paths(config, library_name="L1")
    assert paths == expected
