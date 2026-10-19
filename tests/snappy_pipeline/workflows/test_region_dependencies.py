# -*- coding: utf-8 -*-
"""Tests for the config models of steps that read BED files of regions from ``depends_on``."""

import pydantic
import pytest

from snappy_pipeline.workflows.somatic_purity_ploidy_estimate.model import (
    SomaticPurityPloidyEstimate,
)
from snappy_pipeline.workflows.sv_calling_targeted.model import SvCallingTargeted
from snappy_pipeline.workflows.sv_calling_wgs.model import SvCallingWgs

MELT = {"jar_file": "MELT.jar", "me_refs_path": "me_refs", "me_refs_infix": "1KGP_Hg19"}
SV_DEPENDS_ON = {"alignments": "mapping", "reference": "genome"}


def test_ascat_requires_snp_loci():
    depends_on = {"alignments": "mapping", "reference": "genome"}
    with pytest.raises(pydantic.ValidationError, match="tool=ascat needs depends_on.snp_loci"):
        SomaticPurityPloidyEstimate(tool="ascat", depends_on=depends_on)
    SomaticPurityPloidyEstimate(tool="ascat", depends_on=depends_on | {"snp_loci": "loci"})


@pytest.mark.parametrize("model", [SvCallingTargeted, SvCallingWgs])
def test_melt_requires_gene_regions(model):
    with pytest.raises(pydantic.ValidationError, match="tool=melt needs depends_on.gene_regions"):
        model(tool="melt", melt=MELT, depends_on=SV_DEPENDS_ON)
    model(tool="melt", melt=MELT, depends_on=SV_DEPENDS_ON | {"gene_regions": "genes"})


@pytest.mark.parametrize("model", [SvCallingTargeted, SvCallingWgs])
def test_delly2_excluded_regions_are_optional(model):
    config = model(tool="delly2", delly2={}, depends_on=SV_DEPENDS_ON)
    assert config.depends_on.excluded_regions == ""
