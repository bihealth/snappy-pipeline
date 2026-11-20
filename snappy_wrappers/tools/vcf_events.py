# -*- coding: utf-8 -*-
"""The events of a varlociraptor VCF, and the ones a configuration selects.

A configured entry is an event name, which must exist,
or a shell-style pattern such as ``denovo_*``,
which selects every matching event and may match none:
a pedigree without a child has no de novo event.
"""

from __future__ import annotations

import fnmatch
import gzip
import re

_PROB = re.compile(r"##INFO=<ID=PROB_([A-Za-z0-9_]+),")
_EVENT = re.compile(r'Description="Posterior probability for event (\S+) \(PHRED\)"')
_WILDCARDS = re.compile(r"[*?\[]")


def header_events(path: str) -> list[str]:
    """The events of the ``PROB_*`` fields of a VCF's header,
    named as the scenario names them (from the description;
    the ID in lower case for any other description)."""
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


def select_names(patterns: list[str], available: list[str], what: str = "event") -> list[str]:
    """The names of ``available`` that ``patterns`` select,
    in the order of the patterns and then of ``available``, each once.
    A pattern without wildcards names one of ``available`` and raises ``ValueError`` otherwise;
    a pattern with wildcards (``*``, ``?``, ``[...]``) may match none."""
    selected: list[str] = []
    for pattern in patterns:
        matches = [name for name in available if fnmatch.fnmatchcase(name, pattern)]
        if not matches and not _WILDCARDS.search(pattern):
            raise ValueError(
                f"{what} {pattern!r} is not one of {', '.join(available) or 'none'}; "
                f"a pattern such as {pattern}_* selects every matching {what} and may match none"
            )
        selected += [match for match in matches if match not in selected]
    return selected


def resolve_events(configured: list[str], available: list[str]) -> list[str]:
    """The events of ``available`` that the configured names and patterns select (:func:`select_names`);
    at least one must remain."""
    resolved = select_names(configured, available)
    if not resolved:
        raise ValueError(
            f"none of the events {', '.join(configured)} is in the VCF, which has "
            + ", ".join(available)
        )
    return resolved
