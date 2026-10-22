"""Tests for the tool settings that the ``variant_phasing`` step parts pass to their wrappers."""

from types import SimpleNamespace

from snappy_pipeline.workflows.variant_phasing import (
    PhaseByTransmissionStepPart,
    ReadBackedPhasingStepPart,
)
from snappy_pipeline.workflows.variant_phasing.model import (
    GatkPhaseByTransmission,
    GatkReadBackedPhasing,
)


def _params(part_class, **sections):
    # The params only read the tool's config section, so skip the full step setup.
    part = object.__new__(part_class)
    part.config = SimpleNamespace(**sections)
    return part._get_params_run(SimpleNamespace())


def test_phase_by_transmission_passes_its_settings():
    params = _params(
        PhaseByTransmissionStepPart,
        gatk_phase_by_transmission=GatkPhaseByTransmission(pedigree_validation_type="STRICT"),
    )

    assert params == {
        "de_novo_prior": 1e-8,
        "pedigree_validation_type": "STRICT",
        "father_allele_first": True,
    }


def test_read_backed_phasing_passes_its_settings():
    params = _params(
        ReadBackedPhasingStepPart,
        gatk_read_backed_phasing=GatkReadBackedPhasing(window_size=5000),
    )

    assert params == {
        "phase_quality_threshold": 20.0,
        "window_size": 5000,
        "max_variants_per_window": 100,
    }
