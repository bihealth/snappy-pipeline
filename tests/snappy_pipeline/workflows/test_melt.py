"""Tests for the params of the MELT step part."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.common.melt import MeltStepPart


def _params(action, step="sv_calling_targeted"):
    # The params only read the MELT config section and the step name.
    part = object.__new__(MeltStepPart)
    part.config = SimpleNamespace(
        melt=SimpleNamespace(me_refs_infix="1KGP_Hg19", exclude_contigs=["hs37d5", "NC_007605"])
    )
    part.parent = SimpleNamespace(name=step)
    return getattr(part, f"_get_params_{action}")(SimpleNamespace(me_type="ALU"))


def test_indiv_analysis_gets_exome_for_targeted_data():
    assert _params("indiv_analysis") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "exclude_contigs": ["hs37d5", "NC_007605"],
        "exome": True,
    }


def test_indiv_analysis_gets_no_exome_for_wgs():
    assert "exome" not in _params("indiv_analysis", step="sv_calling_wgs")


@pytest.mark.parametrize("action", ["group_analysis", "genotype", "make_vcf"])
def test_other_commands_do_not_take_exome(action):
    assert _params(action) == {"me_refs_infix": "1KGP_Hg19", "me_type": "ALU"}
