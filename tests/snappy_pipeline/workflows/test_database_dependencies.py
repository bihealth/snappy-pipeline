# -*- coding: utf-8 -*-
"""Tests for the config models of steps that read an index or a database from ``depends_on``."""

import pydantic
import pytest

from snappy_pipeline.workflows.gene_expression_quantification.model import (
    GeneExpressionQuantification,
)
from snappy_pipeline.workflows.repeat_expansion.model import RepeatExpansion
from snappy_pipeline.workflows.somatic_gene_fusion_calling.model import SomaticGeneFusionCalling
from snappy_pipeline.workflows.somatic_wgs_cnv_calling.model import SomaticWgsCnvCalling

FUSION_DEPENDS_ON = {"reads": "data_sets", "reference": "genome", "features": "genes"}


def test_salmon_requires_index():
    depends_on = {"reads": "data_sets", "features": "genes"}
    with pytest.raises(pydantic.ValidationError, match="tool=salmon needs depends_on.index"):
        GeneExpressionQuantification(tool="salmon", salmon={}, depends_on=depends_on)
    GeneExpressionQuantification(
        tool="salmon", salmon={}, depends_on=depends_on | {"index": "salmon_index"}
    )


def test_arriba_requires_index_and_database_is_optional():
    with pytest.raises(pydantic.ValidationError, match="tool=arriba needs depends_on.index"):
        SomaticGeneFusionCalling(tool="arriba", arriba={}, depends_on=FUSION_DEPENDS_ON)
    config = SomaticGeneFusionCalling(
        tool="arriba", arriba={}, depends_on=FUSION_DEPENDS_ON | {"index": "star_index"}
    )
    assert config.depends_on.database == ""


def test_star_fusion_requires_database():
    with pytest.raises(
        pydantic.ValidationError, match="tool=star_fusion needs depends_on.database"
    ):
        SomaticGeneFusionCalling(tool="star_fusion", star_fusion={}, depends_on=FUSION_DEPENDS_ON)
    SomaticGeneFusionCalling(
        tool="star_fusion",
        star_fusion={},
        depends_on=FUSION_DEPENDS_ON | {"database": "star_fusion_database"},
    )


def test_control_freec_requires_database():
    depends_on = {"alignments": "mapping"}
    kwargs = {"tool": "control_freec", "control_freec": {"convert": {}}}
    with pytest.raises(
        pydantic.ValidationError, match="tool=control_freec needs depends_on.database"
    ):
        SomaticWgsCnvCalling(depends_on=depends_on, **kwargs)
    SomaticWgsCnvCalling(depends_on=depends_on | {"database": "control_freec_database"}, **kwargs)


def test_expansionhunter_requires_database():
    depends_on = {"alignments": "mapping", "reference": "genome"}
    with pytest.raises(pydantic.ValidationError, match="database"):
        RepeatExpansion(tool="expansionhunter", depends_on=depends_on)
    RepeatExpansion(tool="expansionhunter", depends_on=depends_on | {"database": "catalog"})
