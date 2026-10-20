# -*- coding: utf-8 -*-
"""Tests for the kit helpers of ``BaseStep``: kits of targeted libraries from the sample sheet,
files from an ``external_data`` task with ``kits``."""

import pandas as pd
import pytest

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.orchestration import load_project
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows import abstract
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, Kits


class _DependsOn(SnappyModel):
    kits: Kits


class _Config(SnappyStepModel):
    depends_on: _DependsOn


class _KitConsumer(BaseStep):
    name = "kit_consumer"
    produces = [DataSignature(DataType.QC)]
    config_model_class = _Config


LIBRARIES = pd.DataFrame(
    {
        "library_name": ["wes1", "wes2", "panel1", "wes_typo", "wes_none", "wgs1", "long1"],
        "library_type": ["WES", "WES", "Panel_seq", "WES", "WES", "WGS", "WGS"],
        "library_kit": ["V6", "V6", "Panel1", "V6 ", None, "TruSeq", "PacBio HiFi"],
    }
)


@pytest.fixture
def step(monkeypatch, mocker, tmp_path):
    monkeypatch.setitem(WORKFLOW_REGISTRY, "kit_consumer", _KitConsumer)
    mocker.patch.object(abstract, "build_library_dataframe", return_value=LIBRARIES)
    kits = {
        "V6": {"targets": "/kits/v6.targets.bed", "baits": "/kits/v6.baits.bed"},
        "Panel1": {"targets": "/kits/p1.bed", "baits": "/kits/p1.bed"},
    }
    config = {
        "data_sets": {},
        "tasks": [
            {
                "step": "external_data",
                "name": "kits",
                "config": {"produces": {"type": "regions", "tags": ["kits"]}, "kits": kits},
            },
            {
                "step": "kit_consumer",
                "name": "consumer",
                "config": {"depends_on": {"kits": "kits"}},
            },
        ],
    }
    project = load_project(config, str(tmp_path))
    # Bypass ``BaseStep.__init__``, which needs a Snakemake workflow and sample sheets.
    step = object.__new__(_KitConsumer)
    step.project = project
    step.task_name = "consumer"
    step.w_config = project.model
    step.config = project.task_configs["consumer"]
    step.depends_on = step.config.depends_on
    step.data_set_infos = step.sheets = []
    step._library_dataframes = {}
    return step


def test_only_targeted_libraries_have_kits(step):
    assert step.library_kits == {"wes1": "V6", "wes2": "V6", "panel1": "Panel1", "wes_typo": "V6 "}


def test_library_kit_paths_are_the_files_of_the_library_kit(step):
    paths = step.library_kit_paths("kits", "panel1")
    assert (paths.targets, paths.baits) == ("/kits/p1.bed", "/kits/p1.bed")


def test_check_names_every_library_without_a_known_kit(step):
    step.check_library_kits("kits", ["wes1", "panel1"])
    with pytest.raises(ValueError, match=r"wes_typo \(V6 \), wes_none \(no kit\), wgs1 \(no kit\)"):
        step.check_library_kits("kits", ["wes1", "wes_typo", "wes_none", "wgs1"])


def test_check_without_require_passes_libraries_without_a_kit(step):
    step.check_library_kits("kits", ["wes1", "wes_none", "wgs1", "long1"], require=False)
    with pytest.raises(ValueError, match="wes_typo"):
        step.check_library_kits("kits", ["wes1", "wes_typo", "wgs1"], require=False)


def test_common_kit(step):
    assert step.common_kit("kits", ["wes1", "wes2"]) == "V6"
    with pytest.raises(ValueError, match="one kit, but they have the kits Panel1, V6"):
        step.common_kit("kits", ["wes1", "panel1"])
