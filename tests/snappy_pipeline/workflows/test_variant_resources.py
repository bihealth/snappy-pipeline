# -*- coding: utf-8 -*-
"""Tests for the ``depends_on`` resources of the Mutect2 config models."""

import pydantic
import pytest

from snappy_pipeline.workflows.panel_of_normals.model import PanelOfNormals
from snappy_pipeline.workflows.variant_calling.model import VariantCalling

DEPENDS_ON = {"alignments": "mapping", "reference": "genome"}


def test_contamination_requires_common_variants():
    with pytest.raises(pydantic.ValidationError, match="depends_on.common_variants is required"):
        VariantCalling(
            tool="mutect2",
            mutect2={"contamination": {"enabled": True}},
            depends_on=DEPENDS_ON,
        )


def test_resources_are_optional_for_mutect2_calling():
    VariantCalling(
        tool="mutect2",
        mutect2={"contamination": {"enabled": False}},
        depends_on=DEPENDS_ON,
    )


def test_contamination_accepts_common_variants():
    VariantCalling(
        tool="mutect2",
        mutect2={"contamination": {"enabled": True}},
        depends_on=DEPENDS_ON | {"common_variants": "common"},
    )


def test_mutect2_panel_requires_germline_resource():
    with pytest.raises(pydantic.ValidationError, match="needs depends_on.germline_resource"):
        PanelOfNormals(tool="mutect2", mutect2={}, depends_on=DEPENDS_ON)
    PanelOfNormals(
        tool="mutect2",
        mutect2={},
        depends_on=DEPENDS_ON | {"germline_resource": "gnomad"},
    )
