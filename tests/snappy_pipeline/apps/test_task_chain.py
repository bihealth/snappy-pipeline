# -*- coding: utf-8 -*-
"""Tests of what the task questions read from the step models."""

from snappy_pipeline.apps.impl.task_chain import (
    TaskEntry,
    candidates,
    entries_of,
    entry_for,
    external_data_settings,
    providers,
    required_keys,
    step_tools,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

CONFIG = {
    "tasks": [
        {
            "name": "genome",
            "step": "external_data",
            "config": {"produces": {"type": "raw", "tags": ["reference", "dna"]}},
        },
        {"name": "kits", "step": "external_data", "config": {"produces": {"type": "regions"}}},
        {"name": "mapping", "step": "ngs_mapping", "config": {"tool": "bwa"}},
        {"name": "calling", "step": "variant_calling", "config": {"tool": "freebayes"}},
        {"step": "variant_calling"},  # no name
    ]
}


def test_entries_of_the_config_carry_tool_and_external_data():
    entries = entries_of(CONFIG)

    assert [entry.name for entry in entries] == ["genome", "kits", "mapping", "calling"]
    assert entries[0].produces == DataSignature(DataType.RAW, frozenset({"reference", "dna"}))
    assert entries[1].produces == DataSignature(DataType.REGIONS)
    assert entries[2] == TaskEntry("mapping", "ngs_mapping", "bwa", None)
    assert entries_of(None) == [] and entries_of({"tasks": None}) == []


def test_step_tools_lists_the_enum_of_a_step():
    assert step_tools("ngs_mapping") == ["bwa", "bwa_mem2", "bwa_mem3", "minimap2", "star"]
    assert step_tools("external_data") == []


def test_required_keys_follow_the_tool():
    bwa = required_keys("ngs_mapping", "bwa")
    assert [key.name for key in bwa] == ["reads", "index", "reference"]
    assert bwa[0].description.startswith("FASTQ source")
    assert bwa[2].contract.signature == DataSignature(DataType.RAW, frozenset({"reference", "dna"}))
    # The features and the exons are required for STAR only
    assert [key.name for key in required_keys("ngs_mapping", "star")] == [
        "reads",
        "index",
        "reference",
        "features",
        "exons",
    ]
    # Keys a setting needs (the igv views of variant_report) are not asked for
    assert [key.name for key in required_keys("variant_report", "datavzrd")] == ["variants"]
    assert [key.name for key in required_keys("variant_filtration", "vembrane")] == ["variants"]
    assert [key.name for key in required_keys("variant_filtration", "dkfz")] == [
        "variants",
        "alignments",
        "reference",
    ]
    assert required_keys("external_data", None) == []


def test_candidates_are_the_fitting_tasks_newest_first():
    entries = entries_of(CONFIG)
    keys = {key.name: key for key in required_keys("variant_calling", "freebayes")}

    assert candidates(keys["alignments"], entries) == ["mapping"]
    assert candidates(keys["reference"], entries) == ["genome"]
    reads = required_keys("ngs_mapping", "bwa")[0]
    assert candidates(reads, entries) == ["genome"]  # the only task producing raw data
    variants = required_keys("variant_filtration", "vembrane")[0]
    assert candidates(variants, entries) == ["calling"]


def test_providers_are_external_data_and_the_steps_producing_the_type():
    keys = {key.name: key for key in required_keys("variant_calling", "freebayes")}

    steps = providers(keys["alignments"])
    assert steps[0] == "external_data"
    assert "ngs_mapping" in steps and "alignment_processing" in steps
    assert "variant_calling" not in steps
    assert "reference_download" in providers(keys["reference"])


def test_external_data_settings_name_the_data_and_the_file_keys():
    keys = {key.name: key for key in required_keys("variant_calling", "freebayes")}

    assert external_data_settings(keys["reference"]) == {
        "produces": {"type": "raw", "tags": ["dna", "reference"]},
        "files": {"fasta": None},
    }
    assert external_data_settings(keys["alignments"]) == {
        "produces": {"type": "alignments", "tags": ["dna"]}
    }
    kits = next(key for key in required_keys("ngs_data_qc", "alfred") if key.name == "kits")
    assert external_data_settings(kits) == {
        "produces": {"type": "regions", "tags": ["kits"]},
        "kits": None,
    }
    entry = entry_for("genome", "external_data", None, external_data_settings(keys["reference"]))
    assert entry.produces == DataSignature(DataType.RAW, frozenset({"dna", "reference"}))
