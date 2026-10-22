"""Tests for the tool settings that the ``variant_calling`` step parts pass to their wrappers."""

from types import SimpleNamespace

from snappy_pipeline.workflows.variant_calling import (
    BcftoolsCallStepPart,
    Gatk3HaplotypeCallerStepPart,
    Gatk3UnifiedGenotyperStepPart,
    Gatk4HaplotypeCallerGvcfStepPart,
    Gatk4HaplotypeCallerJointStepPart,
)
from snappy_pipeline.workflows.variant_calling.model import (
    BcftoolsCall,
    Gatk3Hc,
    Gatk3Ug,
    Gatk4HcGvcf,
    Gatk4HcJoint,
)


def _params(part_class, action="run", **sections):
    # The params only read the tool's config section, so skip the full step setup.
    part = object.__new__(part_class)
    part.config = SimpleNamespace(**{k: v for k, v in sections.items()})
    part.parent = SimpleNamespace(get_upstream_paths=lambda name: SimpleNamespace(fasta="ref.fa"))
    return getattr(part, f"_get_params_{action}")(SimpleNamespace())


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


def test_gatk3_hc_passes_its_settings():
    params = _params(Gatk3HaplotypeCallerStepPart, gatk3_hc=Gatk3Hc(sample_ploidy=3, padding=0))

    assert (params["padding"], params["sample_ploidy"]) == (0, 3)


def test_gatk3_ug_passes_its_settings():
    params = _params(Gatk3UnifiedGenotyperStepPart, gatk3_ug=Gatk3Ug())

    assert (params["padding"], params["sample_ploidy"]) == (1000, 2)


def test_gatk4_hc_joint_passes_its_settings():
    params = _params(Gatk4HaplotypeCallerJointStepPart, gatk4_hc_joint=Gatk4HcJoint(padding=500))

    assert params["padding"] == 500
    assert params["annotation_groups"] == ["StandardAnnotation", "StandardHCAnnotation"]


def test_gatk4_hc_gvcf_passes_the_annotation_groups_of_each_action():
    config = Gatk4HcGvcf(genotype_annotation_groups=["StandardAnnotation"])
    groups = {
        action: _params(Gatk4HaplotypeCallerGvcfStepPart, action, gatk4_hc_gvcf=config)[
            "annotation_groups"
        ]
        for action in ("discover", "combine_gvcfs", "genotype")
    }

    assert groups["genotype"] == ["StandardAnnotation"]
    assert "AlleleSpecificAnnotation" in groups["discover"]
    assert "AlleleSpecificAnnotation" in groups["combine_gvcfs"]
