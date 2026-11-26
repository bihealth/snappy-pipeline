"""The files of libraries: files tables, and the discovery that writes them.

A files table (:class:`FilesTable`) is a TSV with one row per read group of a library (reads:
``library``, ``readgroup``, ``left``, ``right``, ...) or one row per library (other data:
``library`` and the output keys, such as ``bam`` and ``bai``). Relative paths are relative to the
table's directory; lines starting with ``#`` are comments (the first one records the scan that
wrote the table).

:class:`ReadDiscovery` finds the files below search paths, in a folder named after the library.
Search patterns are regular expressions that must match the whole path below that folder; their
``readgroup`` group pairs the mates of one sequencing unit, usually a lane. ``snappy scan-files``
writes the tables with it; ``find_files`` finds other per-library files the same way, one file per
pattern key.
"""

from __future__ import annotations

import csv
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


#: The first column of a files table, and the second one of a reads table
LIBRARY, READGROUP_COLUMN = "library", "readgroup"


@dataclass(frozen=True)
class FilesTable:
    """The files of libraries, read from one or more files tables (:func:`read_files_tables`).

    ``keys`` are the file columns (``left``, ``right``, ... for reads; the output keys otherwise).
    """

    #: The file columns
    keys: tuple[str, ...]
    #: Whether the rows are read groups (``reads``) or one row of files per library
    reads: bool
    #: Library -> its read groups, ordered by the path of the left file
    _read_groups: Mapping[str, tuple[ReadGroup, ...]] = field(default_factory=dict)
    #: Library -> key -> path
    _files: Mapping[str, Mapping[str, str]] = field(default_factory=dict)

    @property
    def libraries(self) -> list[str]:
        """The libraries with rows, in table order."""
        return list(self._read_groups if self.reads else self._files)

    def read_groups(self, library: str) -> list[ReadGroup]:
        """Return the read groups of ``library``; ``[]`` without rows."""
        return list(self._read_groups.get(library, ()))

    def files(self, library: str) -> dict[str, str]:
        """Return key -> path of ``library``; ``{}`` without a row."""
        return dict(self._files.get(library, {}))

    def paths(self) -> list[str]:
        """Every path of the tables, for ``snappy check-files``."""
        if self.reads:
            return [
                p for groups in self._read_groups.values() for g in groups for p in g.paths.values()
            ]
        return [p for files in self._files.values() for p in files.values()]


def read_files_tables(
    paths: Iterable[str], reads: bool, mixed_se_pe: bool = False, pathvars=None
) -> FilesTable:
    """Read the files tables ``paths`` into one :class:`FilesTable`.

    Raises ``ValueError`` (with ``path:line``) for a bad header, an empty required cell, a read
    group or library listed twice, a path listed twice, a library that mixes single-end and
    paired-end read groups, and a library without right files when the table has them and
    ``mixed_se_pe`` is not set.
    """
    rows: list[tuple[str, int, dict[str, str]]] = []
    keys: tuple[str, ...] | None = None
    for path in paths:
        table_keys, table_rows = _read_table(path, reads, pathvars)
        if keys is None:
            keys = table_keys
        elif set(table_keys) != set(keys):
            raise ValueError(
                f"{path}: columns {', '.join(table_keys)} differ from {', '.join(keys)} of the "
                "other tables of the task"
            )
        rows += table_rows
    keys = keys or ()
    seen_paths: dict[str, str] = {}
    for path, line, row in rows:
        for key in keys:
            if (p := row.get(key)) and p in seen_paths:
                raise ValueError(f"{path}:{line}: {p} is also listed at {seen_paths[p]}")
            elif p:
                seen_paths[p] = f"{path}:{line}"
    if reads:
        return FilesTable(keys, True, _read_groups_of(rows, keys, mixed_se_pe))
    files: dict[str, dict[str, str]] = {}
    where: dict[str, str] = {}
    for path, line, row in rows:
        library = row[LIBRARY]
        if library in files:
            raise ValueError(
                f"{path}:{line}: library {library!r} is also listed at {where[library]}"
            )
        files[library] = {key: row[key] for key in keys if row.get(key)}
        where[library] = f"{path}:{line}"
    return FilesTable(keys, False, _files=files)


def files_table_keys(paths: Iterable[str], reads: bool) -> set[str] | None:
    """Return the file columns of the files tables ``paths``, ``None`` if one does not exist yet
    (``snappy scan-files`` writes it; the libraries without rows are reported at setup)."""
    keys: set[str] = set()
    for path in paths:
        if not os.path.exists(path):
            return None
        keys |= set(_read_table(path, reads, None)[0])
    return keys


def _read_table(path: str, reads: bool, pathvars) -> tuple[tuple[str, ...], list]:
    directory = os.path.dirname(os.path.abspath(path))
    try:
        with open(path, newline="") as f:
            lines = list(csv.reader(f, delimiter="\t"))
    except OSError as e:
        raise ValueError(f"Files table {path}: {e.strerror}") from None
    numbered = [
        (n, cells)
        for n, cells in enumerate(lines, start=1)
        if cells and not cells[0].startswith("#")
    ]
    if not numbered or not any(numbered[0][1]):
        raise ValueError(f"{path}:1: the files table is empty")
    header_line, header = numbered[0]
    fixed = [LIBRARY, READGROUP_COLUMN, "left"] if reads else [LIBRARY]
    if header[: len(fixed)] != fixed:
        raise ValueError(f"{path}:{header_line}: the columns must start with {', '.join(fixed)}")
    if len(set(header)) != len(header):
        raise ValueError(f"{path}:{header_line}: a column name is repeated")
    keys = tuple(c for c in header if c not in (LIBRARY, READGROUP_COLUMN))
    if not keys:
        raise ValueError(
            f"{path}:{header_line}: the table needs a file column besides {', '.join(fixed)}"
        )
    required = fixed
    rows = []
    for number, cells in numbered[1:]:
        if len(cells) == 1 and not cells[0].strip():
            continue
        if len(cells) != len(header):
            raise ValueError(f"{path}:{number}: {len(cells)} cells, the header has {len(header)}")
        row = dict(zip(header, (c.strip() for c in cells)))
        for column in required:
            if not row[column]:
                raise ValueError(f"{path}:{number}: {column} is empty")
        for key in keys:
            if row[key]:
                row[key] = _table_path(row[key], directory, pathvars)
        rows.append((path, number, row))
    return keys, rows


def _table_path(value: str, directory: str, pathvars) -> str:
    if pathvars is not None:
        value = pathvars.apply(value)
    return os.path.normpath(os.path.join(directory, value))


def _read_groups_of(rows, keys, mixed_se_pe) -> dict[str, tuple[ReadGroup, ...]]:
    by_library: dict[str, dict[str, tuple[str, int, dict[str, str]]]] = {}
    for path, line, row in rows:
        groups = by_library.setdefault(row[LIBRARY], {})
        if row[READGROUP_COLUMN] in groups:
            seen = groups[row[READGROUP_COLUMN]]
            raise ValueError(
                f"{path}:{line}: read group {row[READGROUP_COLUMN]!r} of {row[LIBRARY]!r} is also "
                f"listed at {seen[0]}:{seen[1]}"
            )
        groups[row[READGROUP_COLUMN]] = (path, line, row)
    result = {}
    for library, groups in by_library.items():
        paths = [{k: r[k] for k in keys if r[k]} for _, _, r in groups.values()]
        if "right" in keys:
            with_right = {"right" in p for p in paths}
            if with_right == {True, False}:
                raise ValueError(f"Library {library!r} mixes single-end and paired-end read groups")
            if with_right == {False} and not mixed_se_pe:
                raise ValueError(
                    f"Library {library!r} has no right files; set mixed_se_pe on the task for "
                    "single-end data"
                )
        common = os.path.commonpath([p for files in paths for p in files.values()])
        if os.path.isfile(common):
            common = os.path.dirname(common)
        read_groups = [
            ReadGroup(name, files, {k: os.path.relpath(p, common) for k, p in files.items()})
            for name, files in zip(groups, paths)
        ]
        result[library] = tuple(sorted(read_groups, key=lambda g: g.left))
    return result
