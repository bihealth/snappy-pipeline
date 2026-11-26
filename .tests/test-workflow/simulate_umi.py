#!/usr/bin/env python3
"""Simulate the reads of the UMI test pipeline, ``pipelines/snappy-umi``.

One tumor library of the IDT duplex design (read structure ``8M+T 8M+T``) is simulated over the kit targets of the pedigree pipeline.
Each molecule is a fragment with two UMIs, ``a`` and ``b``, from the whitelist ``UMIS``.
The strand AB reads ``a`` on read 1 and ``b`` on read 2; the strand BA reads ``b`` on read 1 and ``a`` on read 2.
A molecule has at least one AB pair and, with probability ``1 - P_SINGLE_STRAND``, at least one BA pair, so it can form a duplex consensus.
Planted variants are carried by whole molecules (both strands), at the VAFs of ``VAFS``.
The script writes ``raw/``, ``resources/umis.txt``, ``truth.tsv`` and ``truth_molecules.tsv`` into the pipeline directory.
"""

from __future__ import annotations

import bisect
import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import simulate_pedigree as sp  # noqa: E402
from simulate_pedigree import (  # noqa: E402
    COMPLEMENT,
    CONTIG,
    HERE,
    INSERT_MEAN,
    INSERT_SD,
    READ_LENGTH,
    read_reference,
    read_targets,
    windows,
    with_errors,
)

PIPELINE = HERE / "pipelines" / "snappy-umi"
LIBRARY = "umi001-T1-DNA1-WES1"
SEED = 20261126
UMI_LENGTH = 8
TEMPLATE_LENGTH = READ_LENGTH - UMI_LENGTH
#: Molecules per base of the windows, divided by the fragment length
MOLECULE_DEPTH = 15
MEAN_FAMILY = 3
P_SINGLE_STRAND = 0.3
N_VARIANTS = 20
VAFS = (0.5, 0.2, 0.1, 0.05)


def make_umis(n: int = 32, min_distance: int = 3) -> tuple[str, ...]:
    """Draw ``n`` 8-mers that differ pairwise in at least ``min_distance`` bases."""
    rng = random.Random(f"{SEED}:umis")
    umis: list[str] = []
    while len(umis) < n:
        umi = "".join(rng.choices("ACGT", k=UMI_LENGTH))
        if all(sum(x != y for x, y in zip(umi, other)) >= min_distance for other in umis):
            umis.append(umi)
    return tuple(umis)


UMIS = make_umis()


def plant_variants(reference: str, spans, rng: random.Random) -> list[tuple[int, str, str, float]]:
    """Return (0-based pos, ref, alt, VAF) inside the windows, away from their edges and 12 bases apart."""
    margin = sp.FLANK + READ_LENGTH
    inner = [(s + margin, e - margin) for s, e in spans if e - s > 2 * margin + 40]
    chosen: dict[int, str] = {}
    while len(chosen) < N_VARIANTS:
        start, end = rng.choice(inner)
        pos = rng.randrange(start, end)
        if reference[pos] != "N" and all(abs(pos - other) > 12 for other in chosen):
            chosen[pos] = reference[pos]
    return [
        (pos, ref, "ACGT"[("ACGT".index(ref) + 1) % 4], VAFS[i % len(VAFS)])
        for i, (pos, ref) in enumerate(sorted(chosen.items()))
    ]


def family_size(rng: random.Random) -> int:
    return max(1, int(rng.expovariate(1 / MEAN_FAMILY)))


def covered(pos: int, start: int, length: int) -> bool:
    """Whether a read template of the fragment at ``start`` spans ``pos``."""
    offset = pos - start
    return 0 <= offset < TEMPLATE_LENGTH or length - TEMPLATE_LENGTH <= offset < length


def simulate(reference: str, spans, variants, rng: random.Random, folder: Path) -> list[dict]:
    """Write the paired reads; return one record per molecule with its variant carriers."""
    positions = [v[0] for v in variants]
    quality = "I" * READ_LENGTH
    molecules: list[dict] = []
    pairs = 0
    with (
        gzip.open(folder / f"{LIBRARY}.R1.fastq.gz", "wt") as r1,
        gzip.open(folder / f"{LIBRARY}.R2.fastq.gz", "wt") as r2,
    ):
        for start_w, end_w in spans:
            # molecules = depth * window length / mean fragment length (one fragment covers
            # its length in bases, so this many fragments give MOLECULE_DEPTH per base)
            for _ in range(round(MOLECULE_DEPTH * (end_w - start_w) / INSERT_MEAN)):
                length = max(READ_LENGTH + 10, round(rng.gauss(INSERT_MEAN, INSERT_SD)))
                start = rng.randrange(start_w, end_w - length + 1)
                fragment = list(reference[start : start + length])
                carried, seen = [], []
                lo = bisect.bisect_left(positions, start)
                hi = bisect.bisect_left(positions, start + length)
                for pos, _ref, alt, vaf in variants[lo:hi]:
                    carries = rng.random() < vaf
                    if carries:
                        fragment[pos - start] = alt
                    if covered(pos, start, length):
                        seen.append(pos)
                        if carries:
                            carried.append(pos)
                fragment = "".join(fragment)
                reverse = fragment[::-1].translate(COMPLEMENT)
                a, b = rng.sample(UMIS, 2)
                n_ab = family_size(rng)
                n_ba = 0 if rng.random() < P_SINGLE_STRAND else family_size(rng)
                molecules.append(
                    {
                        "start": start + 1,
                        "end": start + length,
                        "a": a,
                        "b": b,
                        "n_ab": n_ab,
                        "n_ba": n_ba,
                        "seen": seen,
                        "carried": carried,
                    }
                )
                for strand, count in (("AB", n_ab), ("BA", n_ba)):
                    first, second = (a, b) if strand == "AB" else (b, a)
                    forward, backward = (fragment, reverse) if strand == "AB" else (reverse, fragment)
                    for _ in range(count):
                        pairs += 1
                        name = f"{LIBRARY}:{pairs}"
                        t1 = with_errors(forward[:TEMPLATE_LENGTH], rng)
                        t2 = with_errors(backward[:TEMPLATE_LENGTH], rng)
                        r1.write(f"@{name}/1\n{first}{t1}\n+\n{quality}\n")
                        r2.write(f"@{name}/2\n{second}{t2}\n+\n{quality}\n")
    return molecules


def truth_tables(variants, molecules) -> tuple[str, str]:
    """Return ``truth.tsv`` and ``truth_molecules.tsv``; positions are 1-based."""
    truth = ["contig\tpos\tref\talt\tvaf\tmolecules\talt_molecules"]
    for pos, ref, alt, vaf in variants:
        seen = sum(pos in m["seen"] for m in molecules)
        alt_count = sum(pos in m["carried"] for m in molecules)
        truth.append(f"{CONTIG}\t{pos + 1}\t{ref}\t{alt}\t{vaf}\t{seen}\t{alt_count}")
    rows = ["molecule\tstart\tend\tumi_a\tumi_b\tn_ab\tn_ba"]
    for i, m in enumerate(molecules, 1):
        rows.append(f"{i}\t{m['start']}\t{m['end']}\t{m['a']}\t{m['b']}\t{m['n_ab']}\t{m['n_ba']}")
    return "\n".join(truth) + "\n", "\n".join(rows) + "\n"


def main() -> None:
    rng = random.Random(SEED)
    reference = read_reference()
    spans = windows(read_targets())
    variants = plant_variants(reference, spans, rng)
    folder = PIPELINE / "raw" / LIBRARY
    folder.mkdir(parents=True, exist_ok=True)
    (PIPELINE / "resources").mkdir(exist_ok=True)
    (PIPELINE / "resources" / "umis.txt").write_text("".join(f"{u}\n" for u in UMIS))
    molecules = simulate(reference, spans, variants, rng, folder)
    truth, truth_molecules = truth_tables(variants, molecules)
    (PIPELINE / "truth.tsv").write_text(truth)
    (PIPELINE / "truth_molecules.tsv").write_text(truth_molecules)
    pairs = sum(m["n_ab"] + m["n_ba"] for m in molecules)
    duplex = sum(m["n_ab"] >= 1 and m["n_ba"] >= 1 for m in molecules)
    print(f"simulated {len(molecules)} molecules, {pairs} read pairs, {duplex} duplex molecules")


if __name__ == "__main__":
    main()
