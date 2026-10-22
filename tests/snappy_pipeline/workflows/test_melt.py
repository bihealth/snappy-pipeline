"""Tests for the params of the MELT step part."""

from types import SimpleNamespace

from snappy_pipeline.workflows.common.melt import MeltStepPart


def _params(action, step="sv_calling_targeted"):
    # The params only read the MELT config section and the step name.
    part = object.__new__(MeltStepPart)
    part.config = SimpleNamespace(
        melt=SimpleNamespace(
            me_refs_infix="1KGP_Hg19",
            exclude_contigs=["hs37d5", "NC_007605"],
            read_length=150,
            max_no_call_percent=100,
        )
    )
    part.parent = SimpleNamespace(name=step)
    return getattr(part, f"_get_params_{action}")(SimpleNamespace(me_type="ALU"))


def test_indiv_analysis_gets_exome_for_targeted_data():
    assert _params("indiv_analysis") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "read_length": 150,
        "exclude_contigs": ["hs37d5", "NC_007605"],
        "exome": True,
    }


def test_indiv_analysis_gets_no_exome_for_wgs():
    assert "exome" not in _params("indiv_analysis", step="sv_calling_wgs")


def test_genotype_takes_neither_read_length_nor_no_call_percent():
    assert _params("genotype") == {"me_refs_infix": "1KGP_Hg19", "me_type": "ALU"}


def test_group_analysis_gets_the_read_length():
    assert _params("group_analysis") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "read_length": 150,
    }


def test_make_vcf_gets_the_no_call_percent():
    assert _params("make_vcf") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "max_no_call_percent": 100,
    }
