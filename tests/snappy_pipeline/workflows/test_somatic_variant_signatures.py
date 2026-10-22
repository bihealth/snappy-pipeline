# -*- coding: utf-8 -*-
"""Tests for the VCF tabulation params of ``somatic_variant_signatures``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.somatic_variant_signatures import (
    DeconstructSigsStepPart,
    TabulateVariantsStepPart,
)
from snappy_pipeline.workflows.somatic_variant_signatures.model import DeconstructSigs


def test_tabulate_vcf_names_the_tumor_also_without_matched_normal():
    # The params do not depend on the sample sheet, so skip the full step setup.
    part = object.__new__(TabulateVariantsStepPart)

    params = part._get_params_run(SimpleNamespace(tumor_library="T1-DNA1-WES1"))

    assert params == {"tumor_library": "T1-DNA1-WES1"}


def test_deconstruct_sigs_passes_its_options():
    part = object.__new__(DeconstructSigsStepPart)
    part.config = SimpleNamespace(deconstruct_sigs=DeconstructSigs(tri_counts_method="exome"))

    params = part._get_params_run(SimpleNamespace(tumor_library="T1"))

    assert params == {"signatures": "cosmic", "tri_counts_method": "exome"}
