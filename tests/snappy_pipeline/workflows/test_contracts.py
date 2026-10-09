# -*- coding: utf-8 -*-
"""Tests for task contracts: data signatures and upstream resolution in ``BaseStep``."""

import pytest

from snappy_pipeline.orchestration import load_project
from snappy_pipeline.workflow_model import ConfigModel
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.link_in.model import LinkIn
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow
from snappy_pipeline.workflows.variant_filtration.model import (
    ExpectedVariantVcf,
    VariantFiltrationDependsOn,
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


def test_supports_signature_checks_produces():
    assert NgsMappingWorkflow.supports_signature(_signature(DataType.ALIGNMENTS, "dna"))
    assert not NgsMappingWorkflow.supports_signature(_signature(DataType.VARIANTS))
    assert NgsMappingWorkflow.supports_signature(None)
    with pytest.raises(ValueError, match="does not support signature"):
        NgsMappingWorkflow.require_signature(_signature(DataType.VARIANTS))


def test_namespaced_path():
    assert BaseStep.namespaced_path("calling", "output/x.vcf.gz") == "tasks/calling/output/x.vcf.gz"


# Upstream resolution ---------------------------------------------------------------------------


def _w_config(*tasks):
    return ConfigModel(
        static_data_config={"reference": {"path": "/refs/genome.fa"}},
        tasks=[{"step": step, "name": name, "config": config} for step, name, config in tasks],
        data_sets={},
    )


#: mapping -> annotation -> filtration; step configs are not validated by these tests.
FILTRATION_PROJECT = _w_config(
    ("ngs_mapping", "mapping", {}),
    ("variant_annotation", "annotation", {}),
    ("variant_filtration", "filtration", {}),
)


def _bare_step(workflow_cls, task_name, w_config, depends_on, config=None):
    """Return a workflow object with only the attributes that upstream resolution reads.

    ``BaseStep.__init__`` needs sample sheets and a Snakemake workflow; resolution does not.
    """
    step = object.__new__(workflow_cls)
    step.task_name = task_name
    step.w_config = w_config
    step.depends_on = depends_on
    step.config = config
    return step


def _filtration(**depends_on):
    depends_on = VariantFiltrationDependsOn(**{"variant": "annotation", **depends_on})
    return _bare_step(VariantFiltrationWorkflow, "filtration", FILTRATION_PROJECT, depends_on)


def test_get_upstream_paths_namespaces_and_wraps_in_schema():
    paths = _filtration(ngs_mapping="mapping").get_upstream_paths("variant", library_name="L1")

    assert paths == ExpectedVariantVcf(
        vcf="tasks/annotation/output/L1/out/L1.vcf.gz",
        vcf_tbi="tasks/annotation/output/L1/out/L1.vcf.gz.tbi",
    )


def test_get_upstream_paths_keeps_wildcards_without_identifiers():
    paths = _filtration(ngs_mapping="mapping").get_upstream_paths("ngs_mapping")

    assert paths == ExpectedAlignments(
        bam="tasks/mapping/output/{library_name}/out/{library_name}.bam",
        bai="tasks/mapping/output/{library_name}/out/{library_name}.bam.bai",
    )


@pytest.mark.parametrize(
    "depends_on, field, error",
    [
        ({}, "unknown_field", "no depends_on field named 'unknown_field'"),
        ({}, "ngs_mapping", r"depends_on.ngs_mapping is empty or unset"),
        ({"variant": "missing"}, "variant", "upstream task 'missing' .* was not found"),
        # The mapping task produces alignments, not variants.
        ({"variant": "mapping"}, "variant", "does not support signature"),
    ],
)
def test_get_upstream_paths_errors(depends_on, field, error):
    with pytest.raises(ValueError, match=error):
        _filtration(**depends_on).get_upstream_paths(field, library_name="L1")


# get_task_config --------------------------------------------------------------------------------


def _mapping_step(*tasks, link_in=""):
    """Return the workflow object of task "mapping" in a loaded project with ``tasks``."""
    mapping = {"tool": "bwa", "bwa": {"path_index": "/refs/genome"}}
    if link_in:
        mapping["depends_on"] = {"link_in": link_in}
    config = {
        "static_data_config": {"reference": {"path": "/refs/genome.fa"}},
        "tasks": [{"step": step, "name": name, "config": c} for step, name, c in tasks]
        + [{"step": "ngs_mapping", "name": "mapping", "config": mapping}],
        "data_sets": {},
    }
    project = load_project(config, "/projects/p1")
    step = object.__new__(NgsMappingWorkflow)
    step.project = project
    step.task_name = "mapping"
    step.step_name = "ngs_mapping"
    step.config = project.task_configs["mapping"]
    step.depends_on = step.config.depends_on
    return step


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
