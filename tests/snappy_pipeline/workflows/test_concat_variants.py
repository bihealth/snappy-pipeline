"""Tests for ``concat_variants``: the tags of the output, the inputs per entity and the
grouping check."""

from types import SimpleNamespace

import pytest
from pydantic import ValidationError
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.concat_variants import (
    BcftoolsConcatStepPart,
    ConcatVariantsWorkflow,
)
from snappy_pipeline.workflows.concat_variants.model import ConcatVariants

SMALL = DataSignature(
    DataType.VARIANTS, frozenset({"germline", "snv", "indel", "varlociraptor", "filtered"})
)
SV = DataSignature(DataType.VARIANTS, frozenset({"germline", "sv", "varlociraptor", "filtered"}))


def _config(**settings):
    return ConcatVariants(depends_on={"variants": ["fdr", "sv_fdr"]}, tool="bcftools", **settings)


def test_the_output_has_the_common_tags_and_every_variant_kind():
    (produced,) = ConcatVariantsWorkflow.task_produces(_config(), {"variants": (SMALL, SV)})

    assert produced == DataSignature(
        DataType.VARIANTS,
        frozenset({"germline", "varlociraptor", "filtered", "snv", "indel", "sv"}),
    )


def test_a_tag_that_one_input_lacks_is_dropped():
    unfiltered = DataSignature(DataType.VARIANTS, frozenset({"germline", "sv", "varlociraptor"}))

    (produced,) = ConcatVariantsWorkflow.task_produces(_config(), {"variants": (SMALL, unfiltered)})

    assert "filtered" not in produced.tags
    assert {"snv", "indel", "sv"} <= produced.tags


def test_the_inputs_are_the_vcfs_of_every_upstream_task_for_the_entity():
    part = object.__new__(BcftoolsConcatStepPart)
    part.config = _config()
    part.parent = SimpleNamespace(
        upstream_task_names=lambda field: ["fdr", "sv_fdr"],
        get_upstream_paths=lambda field, upstream_task, library_name: SimpleNamespace(
            vcf=f"tasks/{upstream_task}/work/{library_name}/out/{library_name}.vcf.gz",
            vcf_tbi=f"tasks/{upstream_task}/work/{library_name}/out/{library_name}.vcf.gz.tbi",
        ),
    )

    inputs = part._get_input_files_run(Wildcards(fromdict={"library_name": "FAM_P1"}))

    assert inputs["vcf"] == [
        "tasks/fdr/work/FAM_P1/out/FAM_P1.vcf.gz",
        "tasks/sv_fdr/work/FAM_P1/out/FAM_P1.vcf.gz",
    ]
    assert inputs["vcf_tbi"][1] == "tasks/sv_fdr/work/FAM_P1/out/FAM_P1.vcf.gz.tbi"
    assert set(part._get_output_files_run()) == {"vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5"}


def test_the_task_needs_two_inputs_and_a_tool():
    with pytest.raises(ValidationError, match="at least 2 items"):
        ConcatVariants(depends_on={"variants": ["fdr"]}, tool="bcftools")
    with pytest.raises(ValidationError):
        ConcatVariants(depends_on={"variants": ["fdr", "sv_fdr"]})


def test_an_input_with_another_grouping_is_rejected(monkeypatch):
    wf = object.__new__(ConcatVariantsWorkflow)
    wf.config = _config()
    wf.project = SimpleNamespace(group_by={"fdr": "related", "sv_fdr": "library"})
    monkeypatch.setattr(
        ConcatVariantsWorkflow, "effective_group_by", property(lambda self: "related")
    )
    monkeypatch.setattr(
        ConcatVariantsWorkflow, "upstream_group_by", lambda self, task: self.project.group_by[task]
    )
    monkeypatch.setattr(
        ConcatVariantsWorkflow, "upstream_task_names", lambda self, field: ["fdr", "sv_fdr"]
    )

    with pytest.raises(ConfigError, match="'sv_fdr' writes with group_by: library"):
        wf.check_config()
