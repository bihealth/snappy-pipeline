# -*- coding: utf-8 -*-
"""Tests for task contracts: data signatures and upstream resolution in ``BaseStep``."""

import pytest

from snappy_pipeline.orchestration import load_project
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.link_in.model import LinkIn
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_filtration.model import ExpectedVariantVcf
from tests.snappy_pipeline.test_orchestration import (
    WORK_DIR,
    _annotation,
    _calling,
    _config,
    _filtration,
    _mapping,
)


def _signature(data_type, *tags):
    return DataSignature(data_type, frozenset(tags))


@pytest.mark.parametrize(
    "provided, required, expected",
    [
        (_signature(DataType.VARIANTS, "somatic", "snv"), _signature(DataType.VARIANTS), True),
        (_signature(DataType.VARIANTS, "somatic"), _signature(DataType.ALIGNMENTS), False),
        # Plain tags must all be present (AND).
        (
            _signature(DataType.VARIANTS, "somatic", "snv"),
            _signature(DataType.VARIANTS, "somatic", "snv"),
            True,
        ),
        (
            _signature(DataType.VARIANTS, "somatic"),
            _signature(DataType.VARIANTS, "somatic", "snv"),
            False,
        ),
        # A tuple needs at least one of its tags (OR).
        (
            _signature(DataType.VARIANTS, "indel"),
            _signature(DataType.VARIANTS, ("snv", "indel")),
            True,
        ),
        (
            _signature(DataType.VARIANTS, "sv"),
            _signature(DataType.VARIANTS, ("snv", "indel")),
            False,
        ),
        # A "-" prefix forbids a tag (NOT).
        (
            _signature(DataType.VARIANTS, "somatic"),
            _signature(DataType.VARIANTS, "-germline"),
            True,
        ),
        (
            _signature(DataType.VARIANTS, "germline"),
            _signature(DataType.VARIANTS, "-germline"),
            False,
        ),
    ],
)
def test_data_signature_satisfies(provided, required, expected):
    assert provided.satisfies(required) is expected


def test_data_signature_str():
    assert str(_signature(DataType.ALIGNMENTS)) == "alignments"
    assert str(_signature(DataType.VARIANTS, "somatic", ("snv", "indel"))) == (
        "variants [snv|indel, somatic]"
    )


def test_namespaced_path():
    assert BaseStep.namespaced_path("calling", "output/x.vcf.gz") == "tasks/calling/output/x.vcf.gz"


def _step(project, task_name):
    """Return the workflow object of ``task_name`` with only the attributes contract code reads.

    ``BaseStep.__init__`` needs sample sheets and a Snakemake workflow; contract code does not.
    """
    task = project.task(task_name)
    step = object.__new__(WORKFLOW_REGISTRY[task.step])
    step.project = project
    step.task_name = task_name
    step.step_name = task.step
    step.w_config = project.model
    step.config = project.task_configs[task_name]
    step.depends_on = getattr(step.config, "depends_on", None)
    return step


# Upstream resolution ---------------------------------------------------------------------------


def _filtration_step(**depends_on):
    tasks = (_mapping(), _calling(), _annotation(), _filtration(**depends_on))
    return _step(load_project(_config(*tasks), WORK_DIR), "filtration")


def test_get_upstream_paths_namespaces_and_wraps_in_schema():
    paths = _filtration_step().get_upstream_paths("variant", library_name="L1")

    assert paths == ExpectedVariantVcf(
        vcf="tasks/annotation/output/L1/out/L1.vcf.gz",
        vcf_tbi="tasks/annotation/output/L1/out/L1.vcf.gz.tbi",
    )


def test_get_upstream_paths_keeps_wildcards_without_identifiers():
    paths = _filtration_step(ngs_mapping="mapping").get_upstream_paths("ngs_mapping")

    assert paths == ExpectedAlignments(
        bam="tasks/mapping/output/{library_name}/out/{library_name}.bam",
        bai="tasks/mapping/output/{library_name}/out/{library_name}.bam.bai",
    )


@pytest.mark.parametrize(
    "field, error",
    [
        ("unknown_field", "no depends_on field named 'unknown_field'"),
        ("ngs_mapping", r"depends_on.ngs_mapping is empty or unset"),
    ],
)
def test_get_upstream_paths_errors(field, error):
    with pytest.raises(ValueError, match=error):
        _filtration_step().get_upstream_paths(field, library_name="L1")


# get_task_config --------------------------------------------------------------------------------


def _mapping_step(*tasks, link_in=""):
    """Return the workflow object of task "mapping" in a loaded project with ``tasks``."""
    mapping = _mapping(depends_on={"link_in": link_in}) if link_in else _mapping()
    return _step(load_project(_config(*tasks, mapping), WORK_DIR), "mapping")


TRIMMED = ("link_in", "trimmed", {"path": "/data/trimmed"})
RAW = ("link_in", "raw", {"path": "/data/raw"})


def test_get_task_config_returns_own_config():
    step = _mapping_step()
    assert step.get_task_config("ngs_mapping") is step.config
    assert step.get_task_config("mapping") is step.config


def test_get_task_config_follows_depends_on():
    step = _mapping_step(TRIMMED, RAW, link_in="trimmed")
    assert step.get_task_config("link_in") == LinkIn(path="/data/trimmed")


def test_get_task_config_does_not_guess_unset_dependencies():
    # "raw" is the only link_in task, but depends_on.link_in is not set.
    with pytest.raises(ValueError, match="depends_on.link_in is not set"):
        _mapping_step(RAW).get_task_config("link_in")


def test_get_task_config_rejects_unknown_fields():
    with pytest.raises(ValueError, match="depends_on.variant is not set; depends_on fields"):
        _mapping_step().get_task_config("variant")
