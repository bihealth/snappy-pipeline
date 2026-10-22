# -*- coding: utf-8 -*-
"""Tests for the steps and tools that compare each tumor with its matched normal: they leave out
tumors without one, with a warning, and tumor-only tools keep them."""

from types import SimpleNamespace

import pandas as pd
import pytest

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

#: P1 has a tumor and a normal, P2 only a tumor
TABLE = pd.DataFrame(
    {
        "library_name": ["P1-T1-DNA1", "P1-N1-DNA1", "P2-T1-DNA1"],
        "extraction_type": ["dna"] * 3,
        "kind": ["cancer"] * 3,
        "role": ["tumor", "normal", "tumor"],
        "donor_name": ["P1", "P1", "P2"],
        "data_set": ["first"] * 3,
        "is_background": [False] * 3,
    }
)

#: Task name -> step, tool and the tumor libraries the task takes
TASKS = {
    "ascat": ("somatic_purity_ploidy_estimate", "ascat", ["P1-T1-DNA1"]),
    "lohhla": ("somatic_hla_loh_calling", None, ["P1-T1-DNA1"]),
    "cnv_checking": ("somatic_cnv_checking", None, ["P1-T1-DNA1"]),
    "hrd": ("homologous_recombination_deficiency", "scarHRD", ["P1-T1-DNA1"]),
    "sequenza": ("somatic_targeted_seq_cnv_calling", "sequenza", ["P1-T1-DNA1"]),
    "control_freec": ("somatic_wgs_cnv_calling", "control_freec", ["P1-T1-DNA1"]),
    "cnvkit": ("somatic_targeted_seq_cnv_calling", "cnvkit", ["P1-T1-DNA1", "P2-T1-DNA1"]),
    "purecn": ("somatic_targeted_seq_cnv_calling", "purecn", ["P1-T1-DNA1", "P2-T1-DNA1"]),
    "cnvkit_wgs": ("somatic_wgs_cnv_calling", "cnvkit", ["P1-T1-DNA1", "P2-T1-DNA1"]),
    "cnvetti": ("somatic_wgs_cnv_calling", "cnvetti", ["P1-T1-DNA1", "P2-T1-DNA1"]),
}


def _config(step, tool):
    # Built without validation: only the library fields and the tool matter
    model = WORKFLOW_REGISTRY[step].config_model_class
    return model.model_construct(**({"tool": tool} if tool else {}))


def _project():
    return SimpleNamespace(
        task=lambda name: SimpleNamespace(step=TASKS[name][0]),
        task_configs={name: _config(step, tool) for name, (step, tool, _) in TASKS.items()},
        library_tables={False: TABLE},
    )


class _Consumer(BaseStep):
    name = "consumer"
    produces = [DataSignature(DataType.QC)]


def _instance(cls, task_name, config):
    # Bypass ``BaseStep.__init__``; the shared library table is already built
    step = object.__new__(cls)
    step.task_name = task_name
    step.config = config
    step.project = _project()
    step._library_dataframes = {}
    return step


@pytest.mark.parametrize("task", TASKS)
def test_matched_pair_tools_take_only_tumors_with_a_matched_normal(task):
    consumer = _instance(_Consumer, "consumer", None)

    assert consumer.task_libraries(task)["library_name"].tolist() == TASKS[task][2]


@pytest.mark.parametrize("task", TASKS)
def test_tumors_without_a_matched_normal_are_named_in_a_warning(task):
    step, tool, expected = TASKS[task]
    instance = _instance(WORKFLOW_REGISTRY[step], task, _config(step, tool))

    warnings = instance.library_warnings()

    if "P2-T1-DNA1" in expected:
        assert warnings == []
    else:
        assert warnings == [
            f"Task {task!r} leaves out the libraries without a matched normal: P2-T1-DNA1"
        ]
