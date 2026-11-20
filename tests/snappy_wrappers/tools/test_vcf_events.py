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


def test_resolve_events_takes_names_and_patterns():
    available = ["present", "denovo_idx", "denovo_sib1", "recessive", "dominant", "artifact"]

    assert resolve_events(["present", "denovo_*"], available) == [
        "present",
        "denovo_idx",
        "denovo_sib1",
    ]
    assert resolve_events(["denovo_idx", "denovo_*"], available) == ["denovo_idx", "denovo_sib1"]
    assert resolve_events(
        ["somatic_tumor_?igh", "*_low"], ["somatic_tumor_high", "somatic_tumor_low"]
    ) == [
        "somatic_tumor_high",
        "somatic_tumor_low",
    ]
    # A pattern may match nothing: a pair without a child has no de novo event
    assert resolve_events(["present", "denovo_*"], ["present", "artifact"]) == ["present"]
    # A name must exist, so a typo or the old kind fails instead of selecting nothing
    with pytest.raises(ValueError, match="event 'denovo' is not one of present, denovo_idx"):
        resolve_events(["denovo"], available)
    with pytest.raises(ValueError, match="none of the events denovo_\\* is in the VCF"):
        resolve_events(["denovo_*"], ["present", "artifact"])
