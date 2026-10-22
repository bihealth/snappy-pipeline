"""Tests for the tool settings that the ``variant_calling`` step parts pass to their wrappers."""

from types import SimpleNamespace

from snappy_pipeline.workflows.variant_calling import BcftoolsCallStepPart
from snappy_pipeline.workflows.variant_calling.model import BcftoolsCall


def _params(part_class, **sections):
    # The params only read the tool's config section, so skip the full step setup.
    part = object.__new__(part_class)
    part.config = SimpleNamespace(**{k: v for k, v in sections.items()})
    part.parent = SimpleNamespace(get_upstream_paths=lambda name: SimpleNamespace(fasta="ref.fa"))
    return part._get_params_run(SimpleNamespace())


def test_bcftools_call_passes_its_settings():
    params = _params(BcftoolsCallStepPart, bcftools_call=BcftoolsCall())

    assert params["padding"] == 1000
    assert params["annotate"] == "FORMAT/AD,FORMAT/DP"
    assert (params["per_sample_mf"], params["redo_baq"], params["variants_only"]) == (
        True,
        True,
        True,
    )


def test_bcftools_call_passes_changed_settings():
    params = _params(
        BcftoolsCallStepPart,
        bcftools_call=BcftoolsCall(padding=0, redo_baq=False, annotate="FORMAT/DP"),
    )

    assert (params["padding"], params["redo_baq"], params["annotate"]) == (0, False, "FORMAT/DP")
