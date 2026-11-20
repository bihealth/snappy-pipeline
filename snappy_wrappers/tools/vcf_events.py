# -*- coding: utf-8 -*-
"""The events of a varlociraptor VCF, and the ones a configuration means.

A configured name is an event of the VCF, or the kind of the per-member events ``<kind>_<member>``, such as ``denovo`` for ``denovo_idx`` and ``denovo_sib1``.
"""

from __future__ import annotations

import gzip
import re

_PROB = re.compile(r"##INFO=<ID=PROB_([A-Za-z0-9_]+),")
_EVENT = re.compile(r'Description="Posterior probability for event (\S+) \(PHRED\)"')


def header_events(path: str) -> list[str]:
    """The events of the ``PROB_*`` fields of a VCF's header, named as the scenario names them (from the description; the ID in lower case for any other description)."""
    events = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            if match := _PROB.match(line):
                named = _EVENT.search(line)
                events.append(named.group(1) if named else match.group(1).lower())
    return events


def resolve_events(configured: list[str], available: list[str]) -> list[str]:
    """The events of ``available`` that ``configured`` means: a name as itself, else every ``<name>_...`` event (a kind); a name without any event is skipped, since a pedigree without a child has no de novo event."""
    resolved: list[str] = []
    for name in configured:
        matches = [e for e in available if e == name] or [
            e for e in available if e.startswith(f"{name}_")
        ]
        resolved += [match for match in matches if match not in resolved]
    if not resolved:
        raise ValueError(
            f"none of the events {', '.join(configured)} is in the VCF, which has "
            + ", ".join(available)
        )
    return resolved
