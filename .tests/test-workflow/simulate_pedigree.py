#!/usr/bin/env python3
"""Simulate the reads of the pedigree test pipeline, ``pipelines/snappy-germline_pedigree``.

The pedigree exercises every event of the germline scenario on known truth:

- ``pgf`` and ``pgm`` are the paternal grandparents; ``pgm`` is not sequenced and joins the
  scenario as the unsequenced relative that connects ``pgf`` and ``fa``.
- ``mgm`` is the maternal grandmother; the mother's father is unknown, so ``mo`` gets an unobserved
  parent, and the mother's status is unknown.
- ``idx`` and ``sib1`` are the children; ``sib2`` is a third child without sequencing data, whom the
  scenario leaves out.
- Affected are ``pgf``, ``fa`` and ``idx``: a dominant allele on the grandfather's haplotype
  ``pgf.0`` segregates with them, and all three are homozygous for a recessive allele that the
  unaffected members carry at most once.
- ``idx`` and ``sib1`` each carry a de novo variant; ``mo`` carries an allele from her unknown
  father, which ``sib1`` inherits.

Every founder haplotype carries random background variants, inherited without recombination, so
every variant on ``pgf.0`` segregates like the dominant allele.
Reads are simulated over the kit targets of ``KIT_WINDOW`` from the e2e reference, at ``DEPTH``.
The script writes ``raw/`` and ``truth.tsv`` into the pipeline directory, and the sample sheet
and the kit BED, which are committed; ``tests/snappy_pipeline/test_pedigree_simulation.py``
checks that the committed files match.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
PIPELINE = HERE / "pipelines" / "snappy-germline_pedigree"
REFERENCE = HERE / "resources" / "refs" / "subregion.fa"
MEDEXOME = HERE / "pipelines" / "snappy-germline_wes" / "resources" / "Exome-MedExome.chr12.bed"

CONTIG = "12"
#: The kit targets in this window of the contig are simulated (278 targets, 119 kb)
KIT_WINDOW = (6_000_000, 7_000_000)
KIT = "SimExome"
FLANK = 150
DEPTH = 40
READ_LENGTH = 100
INSERT_MEAN, INSERT_SD = 300, 30
ERROR_RATE = 0.002
#: Background variants per base of each founder haplotype
BACKGROUND_DENSITY = 1 / 2500
SEED = 20261120
LIBRARY = "{}-N1-DNA1-WES1"


@dataclass(frozen=True)
class Member:
    name: str
    sex: str
    affected: str
    father: str
    mother: str
    sequenced: bool


#: The sheet rows; ``mof``, the mother's father, is not in the sheet
MEMBERS = (
    Member("pgf", "M", "Y", "0", "0", True),
    Member("pgm", "F", "N", "0", "0", False),
    Member("fa", "M", "Y", "pgf", "pgm", True),
    Member("mgm", "F", "N", "0", "0", True),
    Member("mo", "F", "U", "0", "mgm", True),
    Member("idx", "F", "Y", "fa", "mo", True),
    Member("sib1", "M", "N", "fa", "mo", True),
    Member("sib2", "M", "U", "fa", "mo", False),
)
FOUNDERS = ("pgf", "pgm", "mgm", "mof")
#: Member -> its two haplotypes as founder haplotypes (paternal, maternal); no recombination
HAPLOTYPES = {
    "pgf": ("pgf.0", "pgf.1"),
    "pgm": ("pgm.0", "pgm.1"),
    "mgm": ("mgm.0", "mgm.1"),
    "mof": ("mof.0", "mof.1"),
    "fa": ("pgf.0", "pgm.0"),
    "mo": ("mof.0", "mgm.0"),
    "idx": ("pgf.0", "mgm.0"),
    "sib1": ("pgm.0", "mof.0"),
}
#: The planted variants: label -> the founder haplotypes or the member (de novo) that carry it
PLANTED = {
    "dominant": ("pgf.0",),
    "recessive": ("pgf.0", "pgf.1", "pgm.0", "mgm.0"),
    "denovo_idx": ("idx",),
    "denovo_sib1": ("sib1",),
    "unknown_father": ("mof.0",),
}
AFFECTED = tuple(m.name for m in MEMBERS if m.affected == "Y")
UNAFFECTED_OBSERVED = tuple(m.name for m in MEMBERS if m.affected == "N" and m.sequenced)

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_targets() -> list[tuple[int, int]]:
    """Return the (start, end) of the kit targets in ``KIT_WINDOW``."""
    targets = []
    for line in MEDEXOME.read_text().splitlines():
        contig, start, end = line.split("\t")[:3]
        if contig == CONTIG and KIT_WINDOW[0] <= int(start) < KIT_WINDOW[1]:
            targets.append((int(start), int(end)))
    return targets


def kit_bed() -> str:
    return "".join(f"{CONTIG}\t{start}\t{end}\n" for start, end in read_targets())


def sample_sheet() -> str:
    """The germline sample sheet; the libraries are named by the ``only_secondary_id`` scheme."""
    lines = [
        "[Metadata]",
        "schema\tgermline_variants",
        "schema_version\tv1",
        "title\tsimulated_pedigree",
        "description\tThree generations with planted dominant, recessive and de novo variants",
        "",
        "[Custom Fields]",
        "key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern",
        "familyId\tbioEntity\tfamily\tstring\t0\t0\t0\t0\t0",
        "extractionType\tbioSample\textraction type\tstring\t0\t0\t0\t0\t0",
        "libraryKit\tngsLibrary\texome enrichment kit\tstring\t0\t0\t0\t0\t0",
        "",
        "[Data]",
        "familyId\tpatientName\tfatherName\tmotherName\tsex\tisAffected\textractionType\tlibraryType"
        "\tfolderName\tlibraryKit\thpoTerms",
    ]
    for m in MEMBERS:
        library = LIBRARY.format(m.name) if m.sequenced else "."
        lib_type, kit = ("WES", KIT) if m.sequenced else (".", ".")
        lines.append(
            "\t".join(
                [
                    "PED1",
                    m.name,
                    m.father,
                    m.mother,
                    m.sex,
                    m.affected,
                    "DNA",
                    lib_type,
                    library,
                    kit,
                    ".",
                ]
            )
        )
    return "\n".join(lines) + "\n"


def read_reference() -> str:
    sequence = []
    with open(REFERENCE) as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence).upper()


@dataclass(frozen=True)
class Variant:
    pos: int  # 0-based
    ref: str
    alt: str
    label: str


def plant_variants(reference: str, targets, rng: random.Random) -> dict[str, dict[int, Variant]]:
    """Return haplotype or member -> position -> variant: the planted and the background ones.

    Positions lie inside targets, at least 12 bases apart and away from the target ends.
    """
    inner = [(start + 10, end - 10) for start, end in targets if end - start > 40]
    weights = [end - start for start, end in inner]
    taken: set[int] = set()

    def position() -> int:
        while True:
            start, end = rng.choices(inner, weights)[0]
            pos = rng.randrange(start, end)
            if reference[pos] != "N" and all(abs(pos - other) > 12 for other in taken):
                taken.add(pos)
                return pos

    def snv(pos: int, label: str) -> Variant:
        ref = reference[pos]
        return Variant(pos, ref, rng.choice([b for b in "ACGT" if b != ref]), label)

    def indel(pos: int, label: str) -> Variant:
        length = rng.randint(1, 3)
        if rng.random() < 0.5:
            return Variant(pos, reference[pos : pos + 1 + length], reference[pos], label)
        return Variant(
            pos, reference[pos], reference[pos] + "".join(rng.choices("ACGT", k=length)), label
        )

    variants: dict[str, dict[int, Variant]] = {
        f"{founder}.{i}": {} for founder in FOUNDERS for i in (0, 1)
    }
    variants.update({m.name: {} for m in MEMBERS})
    for label, carriers in PLANTED.items():
        variant = snv(position(), label)
        for carrier in carriers:
            variants[carrier][variant.pos] = variant
    target_bases = sum(end - start for start, end in targets)
    for founder in FOUNDERS:
        for i in (0, 1):
            haplotype = f"{founder}.{i}"
            for _ in range(round(target_bases * BACKGROUND_DENSITY)):
                pos = position()
                variant = indel(pos, haplotype) if rng.random() < 0.1 else snv(pos, haplotype)
                variants[haplotype][variant.pos] = variant
    return variants


def genotypes(variants: dict[str, dict[int, Variant]]) -> dict[int, dict[str, int]]:
    """Return position -> member -> number of alternative alleles, for the members with a genome."""
    result: dict[int, dict[str, int]] = {}
    for member, (paternal, maternal) in HAPLOTYPES.items():
        for haplotype in (paternal, maternal, member):
            for pos in variants.get(haplotype, {}):
                result.setdefault(pos, {})
                result[pos][member] = result[pos].get(member, 0) + 1
    return result


def expected_event(genotype: dict[str, int]) -> str:
    """The event of the germline scenario that should hold the posterior mass of a variant with
    these genotypes: ``denovo`` for a carrier whose parents are reference, ``recessive`` when the
    affected members are homozygous and no observed unaffected member is, ``dominant`` when the
    affected members carry it and no observed unaffected member does, else ``present``."""
    parents = {m.name: (m.father, m.mother) for m in MEMBERS}
    for member, count in genotype.items():
        father, mother = parents.get(member, ("0", "0"))
        if count and father != "0" and mother != "0":
            if not genotype.get(father, 0) and not genotype.get(mother, 0):
                return "denovo"
    affected = [genotype.get(m, 0) for m in AFFECTED]
    unaffected = [genotype.get(m, 0) for m in UNAFFECTED_OBSERVED]
    if all(a == 2 for a in affected) and all(u < 2 for u in unaffected):
        return "recessive"
    if all(a >= 1 for a in affected) and all(u == 0 for u in unaffected):
        return "dominant"
    return "present"


def truth_table(variants, genotype_table) -> str:
    sequenced = [m.name for m in MEMBERS if m.sequenced]
    lines = ["\t".join(["pos", "ref", "alt", "label", "event", *sequenced])]
    by_pos = {v.pos: v for haplotype in variants.values() for v in haplotype.values()}
    for pos in sorted(by_pos):
        variant = by_pos[pos]
        genotype = genotype_table[pos]
        if not any(genotype.get(m) for m in sequenced):
            continue  # on an untransmitted haplotype of an unsequenced founder: no reads carry it
        lines.append(
            "\t".join(
                [
                    str(pos + 1),
                    variant.ref,
                    variant.alt,
                    variant.label,
                    expected_event(genotype),
                    *(str(genotype.get(m, 0)) for m in sequenced),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def windows(targets) -> list[tuple[int, int]]:
    """Return the merged read windows: the targets with ``FLANK`` plus a read length on each side."""
    pad = FLANK + READ_LENGTH
    merged: list[list[int]] = []
    for start, end in sorted(targets):
        if merged and start - pad <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end + pad)
        else:
            merged.append([start - pad, end + pad])
    return [(s, e) for s, e in merged]


def haplotype_sequence(
    reference: str, window: tuple[int, int], variants: dict[int, Variant]
) -> str:
    """The reference in ``window`` with ``variants`` applied."""
    start, end = window
    pieces, cursor = [], start
    for pos in sorted(p for p in variants if start <= p < end):
        variant = variants[pos]
        if pos < cursor or pos + len(variant.ref) > end:
            continue
        pieces.append(reference[cursor:pos])
        pieces.append(variant.alt)
        cursor = pos + len(variant.ref)
    pieces.append(reference[cursor:end])
    return "".join(pieces)


def simulate_reads(member: str, reference, targets, variants, rng: random.Random, out_dir: Path):
    """Write the paired reads of ``member``: fragments of its two haplotypes over the windows."""
    haplotypes = [dict(variants[h]) for h in HAPLOTYPES[member]]
    for h in haplotypes:
        h.update(variants.get(member, {}))  # the de novo variants, on both haplotypes
    library = LIBRARY.format(member)
    folder = out_dir / library
    folder.mkdir(parents=True, exist_ok=True)
    quality = "I" * READ_LENGTH
    with (
        gzip.open(folder / f"{library}.R1.fastq.gz", "wt") as r1,
        gzip.open(folder / f"{library}.R2.fastq.gz", "wt") as r2,
    ):
        count = 0
        for window in windows(targets):
            for haplotype in haplotypes:
                sequence = haplotype_sequence(reference, window, haplotype)
                fragments = round(DEPTH / 2 * len(sequence) / (2 * READ_LENGTH))
                for _ in range(fragments):
                    length = min(
                        len(sequence),
                        max(READ_LENGTH + 10, round(rng.gauss(INSERT_MEAN, INSERT_SD))),
                    )
                    start = rng.randrange(0, len(sequence) - length + 1)
                    fragment = sequence[start : start + length]
                    count += 1
                    name = f"{library}:{count}"
                    r1.write(
                        f"@{name}/1\n{with_errors(fragment[:READ_LENGTH], rng)}\n+\n{quality}\n"
                    )
                    reverse = fragment[::-1].translate(COMPLEMENT)
                    r2.write(
                        f"@{name}/2\n{with_errors(reverse[:READ_LENGTH], rng)}\n+\n{quality}\n"
                    )


def with_errors(read: str, rng: random.Random) -> str:
    bases = list(read)
    for i, base in enumerate(bases):
        if rng.random() < ERROR_RATE:
            bases[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(bases)


def main() -> None:
    targets = read_targets()
    rng = random.Random(SEED)
    reference = read_reference()
    variants = plant_variants(reference, targets, rng)
    PIPELINE.mkdir(parents=True, exist_ok=True)
    (PIPELINE / "resources").mkdir(exist_ok=True)
    (PIPELINE / "resources" / f"{KIT}.chr12.bed").write_text(kit_bed())
    (PIPELINE / "samplesheet.tsv").write_text(sample_sheet())
    (PIPELINE / "truth.tsv").write_text(truth_table(variants, genotypes(variants)))
    for member in MEMBERS:
        if member.sequenced:
            simulate_reads(
                member.name,
                reference,
                targets,
                variants,
                random.Random(f"{SEED}:{member.name}"),
                PIPELINE / "raw",
            )
    print(
        f"simulated {sum(1 for m in MEMBERS if m.sequenced)} libraries over {len(targets)} targets into {PIPELINE}"
    )


if __name__ == "__main__":
    main()
