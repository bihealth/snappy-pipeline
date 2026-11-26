#!/usr/bin/env python3
"""Check the consensus BAMs and metrics of the UMI pipeline against the simulated truth.

``truth_molecules.tsv`` (from ``../../simulate_umi.py``) lists every molecule with its strand family sizes.
The check requires, for the tasks ``umi_duplex`` and ``umi_simplex``:

- the BAM with its ``.bai`` and ``.md5``, coordinate-sorted, every read group with ``SM:<LIBRARY>``;
- the duplex yield metrics: ``ds_duplexes`` at fraction 1.0 within ``TOLERANCE`` of the molecules with both strands;
- the simplex yield metrics, and the row of the library in ``qc_summary.tsv``.

Exit code 1 on the first failure.
Run it from the repository root with ``pixi run python``; it needs only the standard library.
"""

from __future__ import annotations

import csv
import gzip
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIBRARY = "umi001-T1-DNA1-WES1"
TOLERANCE = 0.25


def fail(path: Path, message: str) -> None:
    shown = path.relative_to(HERE) if path.is_relative_to(HERE) else path
    raise SystemExit(f"FAIL {shown}: {message}")


def library_dir(task: str) -> Path:
    return HERE / "tasks" / task / "output" / LIBRARY


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        fail(path, "expected the file to exist")
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bam_header(path: Path) -> str:
    """The header text of a BAM; BGZF is a series of gzip members."""
    with gzip.open(path, "rb") as f:
        if f.read(4) != b"BAM\x01":
            fail(path, "expected the BAM magic")
        (l_text,) = struct.unpack("<i", f.read(4))
        return f.read(l_text).decode()


def check_bam(task: str) -> None:
    bam = library_dir(task) / "out" / f"{LIBRARY}.bam"
    for path in (bam, Path(f"{bam}.bai"), Path(f"{bam}.md5")):
        if not path.exists():
            fail(path, "expected the file to exist")
    lines = bam_header(bam).splitlines()
    hd = [x for x in lines if x.startswith("@HD")]
    if not hd or "SO:coordinate" not in hd[0]:
        fail(bam, f"expected SO:coordinate in the @HD line, got {hd}")
    groups = [x for x in lines if x.startswith("@RG")]
    if not groups:
        fail(bam, "expected at least one @RG line")
    for line in groups:
        if f"SM:{LIBRARY}" not in line.split("\t"):
            fail(bam, f"expected SM:{LIBRARY} in every @RG line, got {line!r}")
    print(f"{task}: {bam.relative_to(HERE)}, {len(groups)} read group(s), coordinate-sorted")


def truth_duplexes() -> int:
    rows = read_table(HERE / "truth_molecules.tsv")
    return sum(int(r["n_ab"]) >= 1 and int(r["n_ba"]) >= 1 for r in rows)


def metrics_path(task: str, kind: str) -> Path:
    name = f"{LIBRARY}.consensus.{kind}_yield_metrics.txt"
    return library_dir(task) / "report" / "umi_qc" / name


def check_duplex_yield() -> None:
    path = metrics_path("umi_duplex", "duplex")
    rows = read_table(path)
    if not rows or float(rows[-1]["fraction"]) != 1.0:
        fail(path, "expected a last row with fraction 1.0")
    found, expected = int(rows[-1]["ds_duplexes"]), truth_duplexes()
    print(f"ds_duplexes: {found}, simulated molecules with both strands: {expected}")
    if abs(found - expected) > TOLERANCE * expected:
        fail(path, f"expected ds_duplexes within {TOLERANCE:.0%} of {expected}, got {found}")


def check_simplex_yield() -> None:
    path = metrics_path("umi_simplex", "simplex")
    rows = read_table(path)
    if not rows:
        fail(path, "expected a header and at least one row")
    print(f"simplex yield, last row: {rows[-1]}")


def check_qc_summary() -> None:
    path = HERE / "tasks" / "umi_duplex" / "output" / "qc_summary" / "out" / "qc_summary.tsv"
    rows = [r for r in read_table(path) if LIBRARY in r.values()]
    if len(rows) != 1:
        fail(path, f"expected one row for {LIBRARY}, found {len(rows)}")
    print({k: v for k, v in rows[0].items() if "duplex" in k.lower()})


def main() -> int:
    for task in ("umi_duplex", "umi_simplex"):
        check_bam(task)
    check_duplex_yield()
    check_simplex_yield()
    check_qc_summary()
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
