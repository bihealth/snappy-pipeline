#!/usr/bin/env python3
"""Check the FDR-controlled calls of the pedigree pipeline against the planted truth.

``truth.tsv`` (from ``../../simulate_pedigree.py``) names every simulated variant with the event
that should hold its posterior mass.
The check reads the VCF of the ``fdr`` task and requires:

- every planted site (the labels without a dot) is called, and its event (``denovo_<member>``,
  ``recessive``, ``dominant`` or ``present``) has the smallest PHRED value of the events, below
  ``MAX_PHRED``;
- at least ``MIN_CALLED`` of all simulated variants are called, and at least ``MIN_CORRECT`` of the
  called ones have their expected event in front;
- no called variant has ``denovo`` or ``recessive`` in front unless the truth says so;
- the calls report of the ``report`` task has its IGV page (``igv/top.html``) below the size
  limit, with igv.js inlined and the planted de novo site of the index among its sites.

Exit code 1 on failure, with the offending sites printed.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
#: The events of the scenario, read from the VCF header: present, denovo_<member>, recessive, dominant
EVENTS: list[str] = []
MAX_PHRED = 10.0
MIN_CALLED = 0.9
MIN_CORRECT = 0.9


IGV_LIMIT_MB = 50


def check_igv_page(truth: list[dict[str, str]]) -> list[str]:
    """The IGV page of the calls report exists, is self-contained and shows the de novo site."""
    pages = sorted((HERE / "tasks" / "report" / "work").glob("*/report/calls/html/igv/top.html"))
    if len(pages) != 1:
        return [f"expected one IGV page of the report task, found {[str(p) for p in pages]}"]
    page = pages[0]
    failures = []
    size = page.stat().st_size
    if size > IGV_LIMIT_MB * 1024 * 1024:
        failures.append(f"{page} has {size} bytes, above {IGV_LIMIT_MB} MB")
    html = page.read_text(encoding="utf-8")
    if "cdn.jsdelivr.net/npm/igv@" in html or "sessionDictionary" not in html:
        failures.append(f"{page}: igv.js is not inlined or the page has no sessions")
    denovo = [row["pos"] for row in truth if row["label"] == "denovo_idx"]
    if denovo and not any(
        f'"{pos}"' in html or f"\t{pos}\t" in html or f">{pos}<" in html for pos in denovo
    ):
        failures.append(f"{page}: the planted de novo site {denovo} is not among its sites")
    print(f"IGV page: {page.relative_to(HERE)}, {size / 1024 / 1024:.1f} MB")
    return failures


def read_truth() -> list[dict[str, str]]:
    lines = (HERE / "truth.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def read_calls() -> dict[tuple[int, str, str], dict[str, float]]:
    """Return (pos, ref, alt) -> event -> PHRED-scaled posterior of the FDR output."""
    vcfs = sorted((HERE / "tasks" / "fdr" / "output").glob("*/out/*.vcf.gz"))
    if len(vcfs) != 1:
        raise SystemExit(f"expected one FDR output VCF, found {[str(v) for v in vcfs]}")
    calls = {}
    with gzip.open(vcfs[0], "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=PROB_"):
                event = line.split("ID=PROB_", 1)[1].split(",", 1)[0].lower()
                if event not in ("artifact", "absent"):
                    EVENTS.append(event)
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            probs = {
                event: float(info.get(f"PROB_{event.upper()}", "inf").split(",")[0])
                for event in EVENTS
            }
            for alt in fields[4].split(","):
                calls[(int(fields[1]), fields[3], alt)] = probs
    return calls


def find_call(row: dict[str, str], calls: dict) -> dict[str, float] | None:
    """The call of a truth variant: the same site, or for an indel a record within 3 bases with
    the same length change (the candidates are left-aligned)."""
    position = int(row["pos"])
    key = (position, row["ref"], row["alt"])
    if key in calls:
        return calls[key]
    shift = len(row["alt"]) - len(row["ref"])
    if shift == 0:
        return None
    for (pos, ref, alt), probs in calls.items():
        if abs(pos - position) <= 3 and len(alt) - len(ref) == shift:
            return probs
    return None


def main() -> int:
    truth, calls = read_truth(), read_calls()
    failures: list[str] = []
    called = correct = 0
    print(f"{len(calls)} calls, {len(truth)} simulated variants")
    print("pos\tlabel\texpected\tfront\t" + "\t".join(EVENTS))
    for row in truth:
        probs = find_call(row, calls)
        planted = "." not in row["label"]
        if probs is None:
            if planted:
                failures.append(f"{row['pos']} {row['label']}: not called")
            continue
        called += 1
        front = min(EVENTS, key=lambda e: probs[e])
        ok = front == row["event"]
        correct += ok
        if planted or not ok:
            values = "\t".join(f"{probs[e]:.1f}" for e in EVENTS)
            print(f"{row['pos']}\t{row['label']}\t{row['event']}\t{front}\t{values}")
        if planted and not (ok and probs[front] < MAX_PHRED):
            failures.append(
                f"{row['pos']} {row['label']}: expected {row['event']}, got {front} at PHRED "
                f"{probs[front]:.1f}"
            )
        if not ok and (front.startswith("denovo") or front == "recessive"):
            failures.append(
                f"{row['pos']} {row['label']}: {front} in front, expected {row['event']}"
            )
    print(f"called {called}/{len(truth)}, expected event in front for {correct}/{called}")
    if called < MIN_CALLED * len(truth):
        failures.append(f"only {called} of {len(truth)} simulated variants called")
    if called and correct < MIN_CORRECT * called:
        failures.append(
            f"only {correct} of {called} called variants have their expected event in front"
        )
    failures += check_igv_page(truth)
    for failure in failures:
        print("FAIL", failure)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
