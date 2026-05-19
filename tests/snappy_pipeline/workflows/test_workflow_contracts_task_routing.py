# -*- coding: utf-8 -*-
"""Contract tests for task routing semantics in BaseStep."""

from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow

from .conftest import patch_module_fs


def _make_base_config(tasks):
    return {
        "static_data_config": {"reference": {"path": "/path/to/ref.fa"}},
        "tasks": tasks,
        "data_sets": {
            "first_batch": {
                "file": "sheet.tsv",
                "search_patterns": [{"left": "*/*/*_R1.fastq.gz", "right": "*/*/*_R2.fastq.gz"}],
                "search_paths": ["/path"],
                "type": "matched_cancer",
                "naming_scheme": "only_secondary_id",
            }
        },
    }


def _mapping_task(name):
    return {
        "step": "ngs_mapping",
        "name": name,
        "config": {"tool": "bwa", "bwa": {"path_index": "/path/to/bwa/index.fa"}},
    }


def _build_workflow(
    dummy_workflow,
    config,
    config_lookup_paths,
    config_paths,
    work_dir,
    cancer_sheet_fake_fs,
    mocker,
    task_name,
):
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    return HlaTypingWorkflow(
        dummy_workflow,
        config,
        config_lookup_paths,
        config_paths,
        work_dir,
        task_name=task_name,
    )


def test_get_task_config_uses_config_level_depends_on(
    dummy_workflow,
    config_lookup_paths,
    config_paths,
    work_dir,
    cancer_sheet_fake_fs,
    mocker,
):
    config = _make_base_config(
        [
            _mapping_task("explicit_map"),
            _mapping_task("typed_map"),
            {"step": "link_in", "name": "preprocessed_fastq", "config": {"path": "/preprocess"}},
            {
                "step": "hla_typing",
                "name": "hla_explicit",
                "depends_on": {"ngs_mapping": "explicit_map", "link_in": "preprocessed_fastq"},
                "config": {
                    "depends_on": {"ngs_mapping": "typed_map"},
                    "tool": "optitype",
                    "optitype": {"max_reads": 5000},
                },
            },
        ]
    )
    step = _build_workflow(
        dummy_workflow,
        config,
        config_lookup_paths,
        config_paths,
        work_dir,
        cancer_sheet_fake_fs,
        mocker,
        task_name="hla_explicit",
    )

    assert str(step.get_task_config("ngs_mapping").tool) == "bwa"
    assert step.get_task_config("ngs_mapping").bwa.path_index == "/path/to/bwa/index.fa"
    assert step.get_preprocessed_path() == "/preprocess"


def test_get_task_config_falls_back_to_typed_depends_on(
    dummy_workflow,
    config_lookup_paths,
    config_paths,
    work_dir,
    cancer_sheet_fake_fs,
    mocker,
):
    config = _make_base_config(
        [
            _mapping_task("typed_map"),
            {
                "step": "hla_typing",
                "name": "hla_typed",
                "config": {
                    "depends_on": {"ngs_mapping": "typed_map"},
                    "tool": "optitype",
                    "optitype": {"max_reads": 5000},
                },
            },
        ]
    )
    step = _build_workflow(
        dummy_workflow,
        config,
        config_lookup_paths,
        config_paths,
        work_dir,
        cancer_sheet_fake_fs,
        mocker,
        task_name="hla_typed",
    )

    assert str(step.get_task_config("ngs_mapping").tool) == "bwa"


def test_get_task_config_falls_back_to_literal_task_name(
    dummy_workflow,
    config_lookup_paths,
    config_paths,
    work_dir,
    cancer_sheet_fake_fs,
    mocker,
):
    config = _make_base_config(
        [
            _mapping_task("literal_map"),
            {
                "step": "hla_typing",
                "name": "hla_literal",
                "config": {"tool": "optitype", "optitype": {"max_reads": 5000}},
            },
        ]
    )
    step = _build_workflow(
        dummy_workflow,
        config,
        config_lookup_paths,
        config_paths,
        work_dir,
        cancer_sheet_fake_fs,
        mocker,
        task_name="hla_literal",
    )

    assert str(step.get_task_config("literal_map").tool) == "bwa"
