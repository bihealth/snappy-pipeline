# -*- coding: utf-8 -*-
"""Tests of the event resolution of the FDR control: the events of a VCF header and the names or kinds a configuration means."""

from pathlib import Path

import pytest

from snappy_wrappers.tools.vcf_events import header_events, resolve_events

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"


def test_header_events_are_the_prob_fields_named_as_the_scenario_names_them(tmp_path):
    assert header_events(str(FIXTURES / "calls.vcf")) == ["denovo", "present", "artifact", "absent"]

    path = tmp_path / "calls.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=PROB_DENOVO_IDX,Number=A,Type=Float,Description="Posterior probability for event denovo_idx (PHRED)">\n'
        '##INFO=<ID=PROB_MYEVENT,Number=A,Type=Float,Description="Posterior probability for event MyEvent (PHRED)">\n'
        '##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Posterior probability for any artifact (PHRED)">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )

    assert header_events(str(path)) == ["denovo_idx", "MyEvent", "artifact"]


def test_resolve_events_takes_names_and_the_members_of_a_kind():
    available = ["present", "denovo_idx", "denovo_sib1", "recessive", "dominant", "artifact"]

    assert resolve_events(["present", "denovo"], available) == [
        "present",
        "denovo_idx",
        "denovo_sib1",
    ]
    assert resolve_events(["denovo_idx", "denovo"], available) == ["denovo_idx", "denovo_sib1"]
    # A kind the scenario lacks is skipped: a pair without a child has no de novo event
    assert resolve_events(["present", "denovo"], ["present", "artifact"]) == ["present"]
    assert resolve_events(["somatic_tumor"], ["somatic_tumor_high", "somatic_tumor_low"]) == [
        "somatic_tumor_high",
        "somatic_tumor_low",
    ]
    with pytest.raises(ValueError, match="none of the events denovo is in the VCF"):
        resolve_events(["denovo"], ["present", "artifact"])
