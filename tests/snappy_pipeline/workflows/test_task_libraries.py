# -*- coding: utf-8 -*-
"""Tests for ``BaseStep.task_libraries``: the libraries of another task, from its config."""

from types import SimpleNamespace

import pandas as pd
import pytest

from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.hla_typing.model import HlaTyping
from snappy_pipeline.workflows.variant_calling.model import VariantCalling

TABLE = pd.DataFrame(
    {
        "library_name": ["P1-T1-DNA1", "P1-N1-DNA1", "P1-T1-RNA1", "F1-N1-DNA1"],
        "extraction_type": ["dna", "dna", "rna", "dna"],
        "kind": ["cancer", "cancer", "cancer", "germline"],
        "role": ["tumor", "normal", "tumor", "index"],
        "donor_name": ["P1", "P1", "P1", "F1"],
        "data_set": ["first"] * 4,
    }
)

#: Task name -> step and config, built without validation: only the library fields matter
TASKS = {
    "optitype": ("hla_typing", HlaTyping.model_construct(tool="optitype")),
    "optitype_dna": (
        "hla_typing",
        HlaTyping.model_construct(tool="optitype", library_selection="extraction_type == 'dna'"),
    ),
    "arcashla": ("hla_typing", HlaTyping.model_construct(tool="arcashla")),
    "mutect2": ("variant_calling", VariantCalling.model_construct(tool="mutect2")),
}


class _Consumer(BaseStep):
    name = "consumer"
    produces = [DataSignature(DataType.QC)]


@pytest.fixture
def step():
    # Bypass ``BaseStep.__init__``; the shared library table is already built
    step = object.__new__(_Consumer)
    step.task_name = "consumer"
    step.project = SimpleNamespace(
        task=lambda name: SimpleNamespace(step=TASKS[name][0]),
        task_configs={name: config for name, (_, config) in TASKS.items()},
        library_tables={False: TABLE},
    )
    return step


@pytest.mark.parametrize(
    "task, expected",
    [
        # The default of each sheet kind, within the extraction types that OptiType supports
        ("optitype", ["P1-T1-DNA1", "F1-N1-DNA1"]),
        ("optitype_dna", ["P1-T1-DNA1", "P1-N1-DNA1", "F1-N1-DNA1"]),
        ("arcashla", ["P1-T1-RNA1"]),
        ("mutect2", ["P1-T1-DNA1"]),
    ],
)
def test_task_libraries_follow_the_config_of_the_task(step, task, expected):
    assert step.task_libraries(task)["library_name"].tolist() == expected
