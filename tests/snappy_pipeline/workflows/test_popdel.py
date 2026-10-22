"""Tests for the PopDel params of ``sv_calling_wgs``."""

from types import SimpleNamespace

from snappy_pipeline.workflows.sv_calling_wgs import PopDelStepPart
from snappy_pipeline.workflows.sv_calling_wgs.model import Popdel


def test_profile_gets_the_regions():
    # The params only read the PopDel config section, so skip the full step setup.
    part = object.__new__(PopDelStepPart)
    part.config = SimpleNamespace(popdel=Popdel())

    regions = part._get_params_profile(SimpleNamespace())["profile_regions"]

    assert len(regions) == 22
    assert regions[0] == "chr1:35000000-36000000"
    assert regions[-1] == "chr22:25000000-26000000"
