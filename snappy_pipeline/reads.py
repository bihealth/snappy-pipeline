"""Central discovery of the FASTQ files of libraries (plans.md R1).

A library's files are in a folder named after it (the sample sheet's folder name) somewhere below
a data set's search paths. Search patterns are regular expressions that must match the whole path
below that folder. Their ``readgroup`` group pairs the mates of one sequencing unit, usually a lane.
``find_files`` finds other per-library files the same way, one file per pattern key.
"""

from __future__ import annotations

import difflib
import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

#: Name of the regex group that identifies the read group of a file
READGROUP = "readgroup"


@dataclass(frozen=True)
class ReadGroup:
    """The files of one read group of a library, by mate (``left``, ``right``, ...)."""

    #: Value of the ``readgroup`` group of the search pattern
    name: str
    #: Mate -> absolute path
    paths: Mapping[str, str]
    #: Mate -> path below the library folder, e.g. ``lane1/x.R1.fastq.gz``
    relpaths: Mapping[str, str] = field(default_factory=dict)

    @property
    def left(self) -> str:
        return self.paths["left"]

    @property
    def right(self) -> str | None:
        return self.paths.get("right")

    def rebased(self, directory: str) -> ReadGroup:
        """Return this read group with each file at its relative path below ``directory``."""
        paths = {mate: os.path.join(directory, rel) for mate, rel in self.relpaths.items()}
        return ReadGroup(self.name, paths, self.relpaths)


def compile_search_pattern(pattern: Mapping[str, str], reads: bool = True) -> dict[str, re.Pattern]:
    """Compile the regexes of one search pattern.

    Read patterns (``reads``) need a ``left`` entry, and each regex a ``readgroup`` group.
    """
    if reads and "left" not in pattern:
        raise ValueError(f"Search pattern {dict(pattern)} has no 'left' entry")
    compiled = {}
    for mate, regex in pattern.items():
        if regex is None:
            continue
        compiled[mate] = re.compile(regex)
        if reads and READGROUP not in compiled[mate].groupindex:
            raise ValueError(f"Search pattern {mate}: {regex!r} has no (?P<{READGROUP}>...) group")
    return compiled


class ReadDiscovery:
    """Finds the FASTQ files of libraries; walks each search path at most once."""

    def __init__(self) -> None:
        self._files: dict[str, tuple[str, ...]] = {}

    def files(self, root: str) -> tuple[str, ...]:
        """Return the paths of all files below ``root``, relative to it and sorted."""
        if root not in self._files:
            found = []
            for dirpath, _, filenames in os.walk(root, followlinks=True):
                rel_dir = os.path.relpath(dirpath, root)
                found += [f if rel_dir == "." else f"{rel_dir}/{f}" for f in filenames]
            self._files[root] = tuple(sorted(found))
        return self._files[root]

    def find(
        self,
        roots: Iterable[str],
        folder_name: str,
        patterns: Iterable[Mapping[str, str]],
        single_end: bool = False,
    ) -> list[ReadGroup]:
        """Return the read groups in the folders named ``folder_name`` below ``roots``.

        A folder of that name may be at any depth. Each file below it is matched against every
        pattern; the ``readgroup`` group pairs the mates. Read groups are ordered by the path of
        their left file below the folder. ``single_end`` accepts read groups without a right mate
        when the patterns have one. Returns an empty list when no file matches; raises
        ``ValueError`` for missing mates, duplicate files and libraries that mix single-end and
        paired-end read groups.
        """
        roots = list(dict.fromkeys(roots))
        compiled = [compile_search_pattern(pattern) for pattern in patterns]
        groups: dict[str, dict[str, tuple[str, str]]] = {}
        for root in roots:
            for rel in self.files(root):
                parts = rel.split("/")
                if folder_name not in parts[:-1]:
                    continue
                inner = "/".join(parts[parts.index(folder_name) + 1 :])
                for pattern in compiled:
                    for mate, regex in pattern.items():
                        if match := regex.fullmatch(inner):
                            group = groups.setdefault(match[READGROUP], {})
                            if mate in group:
                                raise ValueError(
                                    f"Read group {match[READGROUP]!r} of {folder_name!r} has two "
                                    f"{mate} files: {group[mate][0]} and {os.path.join(root, rel)}"
                                )
                            group[mate] = (os.path.join(root, rel), inner)
        if not groups:
            return []
        return self._check(folder_name, groups, compiled, single_end)

    @staticmethod
    def _check(folder_name, groups, compiled, single_end) -> list[ReadGroup]:
        result = []
        for name, mates in groups.items():
            if "left" not in mates:
                raise ValueError(f"Read group {name!r} of {folder_name!r} has no left file")
            paths = {mate: path for mate, (path, _) in mates.items()}
            relpaths = {mate: rel for mate, (_, rel) in mates.items()}
            result.append(ReadGroup(name, paths, relpaths))
        if any("right" in pattern for pattern in compiled):
            with_right = {group.right is not None for group in result}
            if with_right == {True, False}:
                raise ValueError(f"Found mixed single-end and paired-end data for {folder_name!r}")
            if with_right == {False} and not single_end:
                raise ValueError(
                    f"Found no right files for {folder_name!r}; set mixed_se_pe in the data set "
                    "for single-end data"
                )
        return sorted(result, key=lambda group: group.relpaths["left"])

    def find_files(
        self, roots: Iterable[str], folder_name: str, patterns: Iterable[Mapping[str, str]]
    ) -> dict[str, str]:
        """Return output key -> path of the one file per key in the folders ``folder_name``.

        For data other than reads: each key's regex must match the whole path below the folder,
        for exactly one file. Raises ``ValueError`` if no file matches or a key matches twice.
        """
        roots = list(dict.fromkeys(roots))
        compiled = [compile_search_pattern(pattern, reads=False) for pattern in patterns]
        found: dict[str, str] = {}
        for root in roots:
            for rel in self.files(root):
                parts = rel.split("/")
                if folder_name not in parts[:-1]:
                    continue
                inner = "/".join(parts[parts.index(folder_name) + 1 :])
                for pattern in compiled:
                    for key, regex in pattern.items():
                        if regex.fullmatch(inner):
                            if key in found:
                                raise ValueError(
                                    f"Two {key} files for {folder_name!r}: {found[key]} and "
                                    f"{os.path.join(root, rel)}"
                                )
                            found[key] = os.path.join(root, rel)
        if not found:
            raise ValueError(self.missing(roots, folder_name, what="files"))
        return found

    def missing(self, roots: Iterable[str], folder_name: str, what: str = "reads") -> str:
        """Return the error message for a library whose files are not below ``roots``."""
        folders = {d for root in roots for rel in self.files(root) for d in rel.split("/")[:-1]}
        hint = difflib.get_close_matches(folder_name, sorted(folders), n=3)
        message = f"Found no {what} of {folder_name!r} below {', '.join(roots)}"
        return message + (f"; similar folders: {', '.join(hint)}" if hint else "")
