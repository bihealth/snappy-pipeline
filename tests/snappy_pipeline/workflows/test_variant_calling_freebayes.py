"""Tests for freebayes candidate calling in ``variant_calling``: the libraries called together,
the settings passed to the wrapper and the contract of the tool."""

from types import SimpleNamespace

import pandas as pd
from snakemake.iocontainers import Wildcards

from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_calling import FreebayesStepPart, VariantCallingWorkflow
from snappy_pipeline.workflows.variant_calling.model import Freebayes, VariantCalling

#: The library rows of a trio with an RNA library of the child, and of a tumor/normal pair
ROWS = {
    "child-DNA1": {"library_name": "child-DNA1", "matched_normal_lib": ""},
    "father-DNA1": {"library_name": "father-DNA1", "matched_normal_lib": ""},
    "mother-DNA1": {"library_name": "mother-DNA1", "matched_normal_lib": ""},
    "child-RNA1": {"library_name": "child-RNA1", "matched_normal_lib": ""},
    "P1-T1-DNA1": {"library_name": "P1-T1-DNA1", "matched_normal_lib": "P1-N1-DNA1"},
    "P2-T1-DNA1": {"library_name": "P2-T1-DNA1", "matched_normal_lib": ""},
}
PRIMARY_DNA = pd.DataFrame({"library_name": ["child-DNA1", "father-DNA1", "mother-DNA1"]})
RELATED = ["father-DNA1", "child-DNA1", "mother-DNA1", "child-RNA1"]


def _alignments(field, **kwargs):
    if field == "kits":
        return SimpleNamespace(targets=f"kits/{kwargs['library_kit']}.bed")
    if field == "reference":
        return SimpleNamespace(fai="ref.fa.fai")
    name = kwargs["library_name"]
    return SimpleNamespace(bam=f"{name}.bam", bai=f"{name}.bam.bai")


def _part(group_by="library", kits="", **settings):
    part = object.__new__(FreebayesStepPart)
    part.config = SimpleNamespace(
        freebayes=Freebayes(**settings), depends_on=SimpleNamespace(kits=kits)
    )
    part.parent = SimpleNamespace(
        effective_group_by=group_by,
        primary_dna_libraries=lambda: PRIMARY_DNA,
        get_related_libraries=lambda name: RELATED,
        library_row=lambda name: ROWS[name],
        get_upstream_paths=_alignments,
        reference_inputs=lambda: {"reference": "ref.fa"},
        common_kit=lambda field, names: "MedExome",
    )
    return part


def test_members_of_a_pedigree_are_its_primary_dna_libraries_with_the_index_first():
    assert _part("related").members("child-DNA1") == ["child-DNA1", "father-DNA1", "mother-DNA1"]


def test_members_of_a_tumor_are_the_tumor_and_its_matched_normal():
    assert _part().members("P1-T1-DNA1") == ["P1-T1-DNA1", "P1-N1-DNA1"]
    assert _part().members("P2-T1-DNA1") == ["P2-T1-DNA1"]


def test_run_calls_the_members_in_the_window_and_the_kit_targets():
    wildcards = Wildcards(fromdict={"library_name": "P1-T1-DNA1", "scatteritem": "1-of-2"})

    inputs = _part(kits="kits")._get_input_files_run(wildcards)

    assert inputs["bam"] == ["P1-T1-DNA1.bam", "P1-N1-DNA1.bam"]
    assert inputs["region"] == "work/scatter/1-of-2.region.bed"
    assert inputs["targets"] == "kits/MedExome.bed"
    assert "targets" not in _part()._get_input_files_run(wildcards)


def test_run_passes_the_candidate_thresholds():
    params = _part(min_alternate_fraction=0.1)._get_params_run(SimpleNamespace())

    assert params == {"min_alternate_fraction": 0.1, "min_alternate_count": 2}


def test_scatter_windows_are_disjoint_and_skip_the_ignored_contigs():
    params = _part(ignore_contigs=["^chrM$", "_decoy$"])._get_params_scatter(SimpleNamespace())

    assert params == {"ignore_contigs": ["^chrM$", "_decoy$"], "padding": 0}


def test_freebayes_produces_candidates_per_library_without_filter_ids():
    config = VariantCalling.model_construct(tool="freebayes")

    assert VariantCallingWorkflow.task_produces(config, {}) == (
        DataSignature(DataType.VARIANTS, frozenset({"candidates", "snv", "indel"})),
    )
    assert VariantCallingWorkflow.default_group_by(config) is None
    assert VariantCallingWorkflow.filter_ids(config) == frozenset()
    assert VariantCallingWorkflow.task_library_requirement(config) is None


def test_the_germline_callers_keep_their_contract():
    config = VariantCalling.model_construct(tool="gatk4_hc_gvcf")

    assert VariantCallingWorkflow.task_produces(config, {}) == (
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel"})),
    )
    assert VariantCallingWorkflow.default_group_by(config) == "related"
    assert VariantCallingWorkflow.filter_ids(config) == frozenset({"LowQual"})
