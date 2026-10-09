# -*- coding: utf-8 -*-
"""Tests for loading a project and selecting targets in ``snappy_pipeline.orchestration``."""

import logging

import pydantic
import pytest

from snappy_pipeline import orchestration
from snappy_pipeline.orchestration import Project, load_project, select_target_tasks
from snappy_pipeline.workflow_model import ConfigModel
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

WORK_DIR = "/projects/p1"


def _config(*tasks):
    return {
        "static_data_config": {"reference": {"path": "/refs/genome.fa"}},
        "tasks": [{"step": step, "name": name, "config": config} for step, name, config in tasks],
        "data_sets": {},
    }


def _mapping(name="mapping", **config):
    config = {"depends_on": {"reads": "data_sets"}, **config}
    return ("ngs_mapping", name, {"tool": "bwa", "bwa": {"path_index": "/refs/genome"}, **config})


def _calling(name="calling", mapping="mapping", tool="mutect2"):
    config = {"depends_on": {"alignments": mapping}, "tool": tool}
    config[tool] = {"contamination": {}} if tool == "mutect2" else {}
    return ("variant_calling", name, config)


def _annotation(name="annotation", variants="calling"):
    return (
        "variant_annotation",
        name,
        {"depends_on": {"variants": variants}, "tool": "vep", "vep": {}},
    )


def _filtration(name="filtration", variants="annotation", **depends_on):
    config = {"depends_on": {"variants": variants, **depends_on}, "tool": "bcftools"}
    config["bcftools"] = {"include": "QUAL > 10"}
    return ("variant_filtration", name, config)


def _tmb(name="tmb", variants="filtration"):
    config = {"depends_on": {"variants": variants}}
    return ("tumor_mutational_burden", name, config | {"target_regions": "/refs/regions.bed"})


# load_project ------------------------------------------------------------------------------------


def test_load_project_resolves_dependencies_and_orders_tasks():
    # Listed out of order on purpose.
    project = load_project(_config(_filtration(), _annotation(), _calling(), _mapping()), WORK_DIR)

    assert [task.name for task in project.tasks] == [
        "mapping",
        "calling",
        "annotation",
        "filtration",
    ]
    assert project.dependencies["filtration"] == {"variants": "annotation"}
    assert project.dependencies["mapping"] == {}
    assert project.task_configs["calling"].tool == "mutect2"
    assert project.lookup_paths == (WORK_DIR, "/projects")
    assert project.config_paths == (f"{WORK_DIR}/config.yaml",)


def test_load_project_rejects_unknown_step():
    with pytest.raises(ValueError, match="unknown step 'variant_caling'"):
        load_project(_config(("variant_caling", "calling", {})), WORK_DIR)


def test_load_project_rejects_duplicate_task_names():
    with pytest.raises(ValueError, match="duplicates: mapping"):
        load_project(_config(_mapping(), _mapping()), WORK_DIR)


def test_load_project_rejects_dependency_on_missing_task():
    with pytest.raises(
        ValueError, match="depends_on.alignments is 'mapping2', which is not a task"
    ):
        load_project(_config(_mapping(), _calling(mapping="mapping2")), WORK_DIR)


def test_load_project_rejects_dependency_cycles():
    tasks = (_annotation(variants="filtration"), _filtration(variants="annotation"))
    with pytest.raises(ValueError, match="Dependency cycle between tasks"):
        load_project(_config(*tasks), WORK_DIR)


def test_load_project_uses_defaults_for_keys_without_value(caplog):
    with caplog.at_level(logging.INFO, logger="snappy_pipeline.orchestration"):
        mapping = _mapping(bwa={"path_index": "/refs/genome", "mask_duplicates": None})
        project = load_project(_config(mapping), WORK_DIR)

    assert project.task_configs["mapping"].bwa.mask_duplicates is True
    assert "tasks[0].config.bwa.mask_duplicates has no value" in caplog.text


def test_load_project_validates_step_configs():
    with pytest.raises(pydantic.ValidationError, match="contamination"):
        load_project(
            _config(("variant_calling", "calling", {"tool": "mutect2", "mutect2": {}})), WORK_DIR
        )


def test_load_project_passes_variant_tags_through_annotation_and_filtration():
    tasks = (_mapping(), _calling(), _annotation(), _filtration(), _tmb())
    project = load_project(_config(*tasks), WORK_DIR)

    somatic = frozenset({"somatic", "snv", "indel"})
    assert project.signatures["mapping"] == (
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
    )
    assert project.signatures["filtration"] == (
        DataSignature(DataType.VARIANTS, somatic | {"annotated", "filtered"}),
    )


def test_load_project_rejects_germline_variants_for_tmb():
    tasks = (_mapping(), _calling(tool="gatk4_hc_gvcf"), _annotation(), _filtration(), _tmb())
    with pytest.raises(
        ValueError,
        match=r"Task 'tmb': depends_on.variants requires variants \[somatic\], "
        r"but task 'filtration' produces variants \[annotated, filtered, germline, indel, snv\]",
    ):
        load_project(_config(*tasks), WORK_DIR)


def test_load_project_rejects_rna_alignments_for_variant_calling(tmp_path):
    for index_file in ("Genome", "SA", "SAindex"):
        (tmp_path / index_file).touch()
    star = _mapping("star", tool="star", star={"path_index": str(tmp_path)})
    with pytest.raises(
        ValueError,
        match=r"Task 'calling': depends_on.alignments requires alignments \[dna\], "
        r"but task 'star' produces alignments \[rna\]",
    ):
        load_project(_config(star, _calling(mapping="star")), WORK_DIR)


def test_load_project_rejects_dna_alignments_for_expression_quantification():
    strandedness = {"tool": "strandedness", "strandedness": {"path_exon_bed": "/refs/exons.bed"}}
    expression = ("gene_expression_quantification", "expression", strandedness)
    expression[2]["depends_on"] = {"alignments": "mapping"}
    with pytest.raises(
        ValueError,
        match=r"Task 'expression': depends_on.alignments requires alignments \[rna\], "
        r"but task 'mapping' produces alignments \[dna\]",
    ):
        load_project(_config(_mapping(), expression), WORK_DIR)


def test_load_project_requires_depends_on_keys():
    # No default task name such as "ngs_mapping" fills in a missing key.
    calling = ("variant_calling", "calling", {"tool": "mutect2", "mutect2": {"contamination": {}}})
    with pytest.raises(pydantic.ValidationError, match="depends_on\n  Field required"):
        load_project(_config(_mapping(), calling), WORK_DIR)
    with pytest.raises(pydantic.ValidationError, match="depends_on.reads\n  Field required"):
        load_project(_config(_mapping(depends_on={})), WORK_DIR)


def test_load_project_accepts_reads_from_data_sets():
    project = load_project(_config(_mapping(depends_on={"reads": "data_sets"})), WORK_DIR)
    assert project.dependencies["mapping"] == {}


def test_load_project_reserves_data_sets_for_reads():
    with pytest.raises(
        ValueError, match="depends_on.alignments is 'data_sets', which is not a task"
    ):
        load_project(_config(_mapping(), _calling(mapping="data_sets")), WORK_DIR)
    with pytest.raises(ValueError, match="Task name 'data_sets' is reserved"):
        load_project(_config(_mapping(name="data_sets")), WORK_DIR)


# create_task_instances ---------------------------------------------------------------------------


def test_create_task_instances_builds_each_task_once(monkeypatch):
    project = load_project(_config(_calling(), _mapping()), WORK_DIR)
    built = []

    class Stub:
        def __init__(self, workflow, project, task_name):
            built.append(task_name)

    for step in ("ngs_mapping", "variant_calling"):
        monkeypatch.setitem(WORKFLOW_REGISTRY, step, Stub)
    instances = orchestration.create_task_instances(object(), project)

    assert built == ["mapping", "calling"]
    assert orchestration.task_instance("calling") is instances["calling"]


# select_target_tasks -----------------------------------------------------------------------------


def _project(**dependencies):
    """Return a project with the given task -> {field: upstream} dependencies."""
    tasks = [{"step": "dummy", "name": name, "config": {}} for name in dependencies]
    model = ConfigModel(
        static_data_config={"reference": {"path": "/refs/genome.fa"}}, tasks=tasks, data_sets={}
    )
    return Project(
        config={},
        model=model,
        work_dir=WORK_DIR,
        lookup_paths=(WORK_DIR,),
        config_paths=(),
        tasks=tuple(model.tasks),
        task_configs={},
        dependencies=dependencies,
        signatures={},
    )


#: mapping -> calling -> annotation -> filtration, plus a QC task on the mapping.
CHAIN = _project(
    mapping={},
    calling={"alignments": "mapping"},
    annotation={"variants": "calling"},
    filtration={"variants": "annotation", "alignments": "mapping"},
    qc={"alignments": "mapping"},
)


def test_default_targets_leaf_tasks_in_config_order():
    assert select_target_tasks(CHAIN) == ["filtration", "qc"]


def test_single_task():
    assert select_target_tasks(CHAIN, target_task="calling") == ["calling"]


def test_unknown_single_task_raises():
    with pytest.raises(ValueError, match="Unknown task 'calls'"):
        select_target_tasks(CHAIN, target_task="calls")


def test_all_tasks():
    assert select_target_tasks(CHAIN, all_tasks=True) == [
        "mapping",
        "calling",
        "annotation",
        "filtration",
        "qc",
    ]
