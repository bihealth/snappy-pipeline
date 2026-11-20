"""The simulated pedigree of the e2e pipeline ``snappy-germline_pedigree``: the committed sample
sheet and kit BED match what ``simulate_pedigree.py`` writes, and the planted variants land in the
events they are meant for."""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / ".tests" / "test-workflow" / "simulate_pedigree.py"
PIPELINE = ROOT / ".tests" / "test-workflow" / "pipelines" / "snappy-germline_pedigree"

spec = importlib.util.spec_from_file_location("simulate_pedigree", SCRIPT)
assert spec is not None and spec.loader is not None
sim = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = sim  # dataclasses resolve the annotations through sys.modules
spec.loader.exec_module(sim)


def test_the_committed_sheet_and_kit_match_the_simulation():
    assert (PIPELINE / "samplesheet.tsv").read_text() == sim.sample_sheet()
    assert (PIPELINE / "resources" / "SimExome.chr12.bed").read_text() == sim.kit_bed()


def test_the_planted_variants_have_their_events():
    variants = {name: {} for name in sim.HAPLOTYPES}
    variants.update({f"{founder}.{i}": {} for founder in sim.FOUNDERS for i in (0, 1)})
    for pos, (label, carriers) in enumerate(sim.PLANTED.items()):
        for carrier in carriers:
            variants[carrier][pos] = sim.Variant(pos, "A", "C", label)
    events = {
        label: sim.expected_event(sim.genotypes(variants)[pos])
        for pos, label in enumerate(sim.PLANTED)
    }
    assert events == {
        "dominant": "dominant",
        "recessive": "recessive",
        "denovo_idx": "denovo_idx",
        "denovo_sib1": "denovo_sib1",
        "unknown_father": "present",
    }


def test_the_pedigree_has_the_constellations_it_claims():
    by_name = {m.name: m for m in sim.MEMBERS}
    assert not by_name["pgm"].sequenced and by_name["fa"].mother == "pgm"  # connecting relative
    assert by_name["mo"].father == "0" and by_name["mo"].mother == "mgm"  # one parent unknown
    assert by_name["mo"].affected == "U" and by_name["sib2"].affected == "U"
    assert not by_name["sib2"].sequenced  # left out
    assert {m.sex for m in sim.MEMBERS} == {"M", "F"}
