# -*- coding: utf-8 -*-
"""Tests for task contracts: data signatures and upstream resolution in ``BaseStep``."""

import pytest

from snappy_pipeline.orchestration import load_project
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, select_signature
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


CALLS = _signature(DataType.VARIANTS, "somatic", "snv", "indel")
PASS_CALLS = CALLS.with_tags("filtered")


@pytest.mark.parametrize(
    "required, expected",
    [
        # PASS calls are the default when both satisfy the requirement.
        (_signature(DataType.VARIANTS, "somatic"), PASS_CALLS),
        (None, PASS_CALLS),
        (_signature(DataType.VARIANTS, "-filtered"), CALLS),
        (_signature(DataType.VARIANTS, "germline"), None),
    ],
)
def test_select_signature_prefers_filtered_calls(required, expected):
    assert select_signature((CALLS, PASS_CALLS), required) == expected


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
    paths = _filtration_step().get_upstream_paths("variants", library_name="L1")

    assert paths == ExpectedVariantVcf(
        vcf="tasks/annotation/output/L1/out/L1.vcf.gz",
        vcf_tbi="tasks/annotation/output/L1/out/L1.vcf.gz.tbi",
    )


def test_get_upstream_paths_keeps_wildcards_without_identifiers():
    paths = _filtration_step(alignments="mapping").get_upstream_paths("alignments")

    assert paths == ExpectedAlignments(
        bam="tasks/mapping/output/{library_name}/out/{library_name}.bam",
        bai="tasks/mapping/output/{library_name}/out/{library_name}.bam.bai",
    )


@pytest.mark.parametrize(
    "field, error",
    [
        ("unknown_field", "no depends_on field named 'unknown_field'"),
        ("alignments", r"depends_on.alignments is empty or unset"),
    ],
)
def test_get_upstream_paths_errors(field, error):
    with pytest.raises(ValueError, match=error):
        _filtration_step().get_upstream_paths(field, library_name="L1")


# get_task_config --------------------------------------------------------------------------------


def _mapping_step(*tasks, reads="data_sets"):
    """Return the workflow object of task "mapping" in a loaded project with ``tasks``."""
    mapping = _mapping(depends_on={"reads": reads})
    return _step(load_project(_config(*tasks, mapping), WORK_DIR), "mapping")


READS_PATTERN = {
    "left": r"(?P<readgroup>.+)\.R1\.fastq\.gz",
    "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz",
}


def _external_reads(name, path):
    config = {"produces": {"type": "raw"}, "search_paths": [path]}
    return ("external_data", name, config | {"search_patterns": [READS_PATTERN]})


TRIMMED = _external_reads("trimmed", "/data/trimmed")
RAW = _external_reads("raw", "/data/raw")


def test_get_task_config_returns_own_config():
    step = _mapping_step()
    assert step.get_task_config("ngs_mapping") is step.config
    assert step.get_task_config("mapping") is step.config


def test_get_task_config_follows_depends_on():
    step = _mapping_step(TRIMMED, RAW, reads="trimmed")
    assert step.get_task_config("reads").search_paths == ["/data/trimmed"]


def test_get_task_config_does_not_guess_unset_dependencies():
    # "raw" is the only external_data task, but depends_on.reads is data_sets.
    with pytest.raises(ValueError, match="depends_on.reads names no task"):
        _mapping_step(RAW).get_task_config("reads")


def test_get_task_config_rejects_unknown_fields():
    with pytest.raises(ValueError, match="depends_on.variant names no task; depends_on fields"):
        _mapping_step().get_task_config("variant")


# Filtered and unfiltered calls ------------------------------------------------------------------


def test_mutect2_calling_provides_pass_calls_by_default_and_all_calls_on_request():
    project = load_project(_config(_mapping(), _calling(), _annotation()), WORK_DIR)
    assert project.signatures["calling"] == (CALLS, PASS_CALLS)
    # Annotation reads the PASS calls, so its output carries the tag on.
    assert project.signatures["annotation"] == (PASS_CALLS.with_tags("annotated"),)

    annotation = _step(project, "annotation")
    default = annotation.get_upstream_paths("variants", library_name="T1")
    unfiltered = annotation.get_upstream_paths(
        "variants", signature=_signature(DataType.VARIANTS, "-filtered"), library_name="T1"
    )
    assert default.vcf == "tasks/calling/output/T1/out/T1.vcf.gz"
    assert unfiltered.vcf == "tasks/calling/output/T1/out/T1.full.vcf.gz"


def test_unfiltered_calls_are_not_available_after_annotation():
    project = load_project(_config(_mapping(), _calling(), _annotation(), _filtration()), WORK_DIR)
    with pytest.raises(ValueError, match=r"produces no variants \[-filtered\]"):
        _step(project, "filtration").get_upstream_paths(
            "variants", signature=_signature(DataType.VARIANTS, "-filtered"), library_name="T1"
        )


# External data ------------------------------------------------------------------------------------


def test_external_files_are_used_where_they_are(tmp_path):
    for name in ("T1.bam", "T1.bam.bai"):
        (tmp_path / "T1").mkdir(exist_ok=True)
        (tmp_path / "T1" / name).touch()
    bams = (
        "external_data",
        "bams",
        {
            "produces": {"type": "alignments", "tags": ["dna"]},
            "search_paths": [str(tmp_path)],
            "search_patterns": [{"bam": r".+\.bam", "bai": r".+\.bam\.bai"}],
        },
    )
    project = load_project(_config(bams, _calling(mapping="bams")), WORK_DIR)

    alignments = _step(project, "calling").get_upstream_paths("alignments", library_name="T1")
    assert alignments.bam == str(tmp_path / "T1" / "T1.bam")
    assert alignments.bai == str(tmp_path / "T1" / "T1.bam.bai")


# Mapper index ---------------------------------------------------------------------------------


def test_mapping_reads_the_index_files_of_its_index_task():
    step = _step(load_project(_config(_mapping()), WORK_DIR), "mapping")
    assert step.get_index_files("bwa") == [
        f"/refs/genome{ext}" for ext in (".amb", ".ann", ".bwt", ".pac", ".sa")
    ]


def test_mapping_rejects_the_index_of_another_tool():
    step = _step(load_project(_config(_mapping()), WORK_DIR), "mapping")
    with pytest.raises(ValueError, match="bwa_index"):
        step.get_index_path("star")
