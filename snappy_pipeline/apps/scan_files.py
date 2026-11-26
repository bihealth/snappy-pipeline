"""``snappy scan-files``: write the files tables of ``external_data`` tasks from the folders.

For every library of the sample sheets, the files are searched in a folder named like the library
(the sheet's ``folder_name``, else the library name) below the search paths, with the search
patterns, as :class:`~snappy_pipeline.reads.ReadDiscovery` does. The rows go to the first table of
the task's ``library_files``; paths below the table's directory are written relative to it, others
absolute. The first line of the table records the command that wrote it, so that a bare
``snappy scan-files`` rescans the same way; without that line, reads are searched below ``raw``
with the default patterns, and other data needs ``--pattern key=regex``.
"""

from __future__ import annotations

import os
import shlex
from collections.abc import Iterable
from dataclasses import dataclass

from snappy_pipeline.library_table import build_library_dataframe
from snappy_pipeline.orchestration import Project
from snappy_pipeline.reads import (
    LIBRARY,
    READGROUP_COLUMN,
    ReadDiscovery,
    compile_search_pattern,
)

#: Where reads are searched without other arguments
DEFAULT_SEARCH_PATHS = ("raw",)
#: The patterns of reads without other arguments
DEFAULT_PATTERNS = (
    {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"},
)
#: The first line of a table written by the command
COMMENT = "# snappy scan-files"


@dataclass(frozen=True)
class Scan:
    """What to search: the search paths (absolute) and the patterns, by key."""

    search_paths: tuple[str, ...]
    patterns: tuple[dict[str, str], ...]

    def command(self, task_name: str, table: str) -> str:
        """The command line that repeats this scan, for the first line of ``table``; search
        paths below the table's directory are relative to it, like the paths of the rows."""
        directory = os.path.dirname(os.path.abspath(table))
        args = ["--task", task_name]
        for path in self.search_paths:
            relative = path.startswith(directory + os.sep)
            args += ["--search-path", os.path.relpath(path, directory) if relative else path]
        for pattern in self.patterns:
            args += [f"--pattern={key}={regex}" for key, regex in pattern.items()]
        return COMMENT + " " + shlex.join(args)


def parse_patterns(values: Iterable[str], reads: bool) -> tuple[dict[str, str], ...]:
    """Turn ``key=regex`` arguments into one pattern; raises ``ValueError`` for a bad one."""
    pattern: dict[str, str] = {}
    for value in values:
        key, sep, regex = value.partition("=")
        if not sep or not key or not regex:
            raise ValueError(f"--pattern takes key=regex, not {value!r}")
        pattern[key] = regex
    if not pattern:
        return ()
    compile_search_pattern(pattern, reads=reads)
    return (pattern,)


def recorded_scan(table: str) -> Scan | None:
    """The scan that the first line of ``table`` records, ``None`` without one."""
    try:
        with open(table) as f:
            first = f.readline()
    except OSError:
        return None
    if not first.startswith(COMMENT):
        return None
    args = shlex.split(first[len(COMMENT) :])
    directory = os.path.dirname(os.path.abspath(table))
    search_paths, patterns, pattern = [], [], {}
    for arg, value in zip(args, args[1:] + [""]):
        if arg == "--search-path":
            search_paths.append(os.path.normpath(os.path.join(directory, value)))
        elif arg.startswith("--pattern="):
            key, _, regex = arg[len("--pattern=") :].partition("=")
            pattern[key] = regex
    if pattern:
        patterns.append(pattern)
    return Scan(tuple(search_paths), tuple(patterns))


def scan_tasks(project: Project, task_names: Iterable[str] = ()) -> list[str]:
    """The ``external_data`` tasks with files tables, or those of ``task_names``; raises
    ``ValueError`` for a name without tables."""
    names = [
        name
        for name, config in project.task_configs.items()
        if getattr(config, "library_files", None)
    ]
    if task_names:
        for name in task_names:
            if name not in names:
                raise ValueError(f"Task {name!r} is no external_data task with library_files")
        return list(task_names)
    return names


def scan_of(project: Project, task_name: str, search_paths=(), patterns=()) -> Scan | None:
    """The scan of a task: the arguments given, else the one its table records, else the defaults
    for reads; ``None`` for other data without arguments or a recorded scan."""
    config = project.task_configs[task_name]
    reads = config.patterns_find_reads()
    recorded = recorded_scan(config.library_files[0])
    if search_paths or patterns:
        given = tuple(os.path.join(project.project_dir, p) for p in search_paths)
        paths = given or (recorded.search_paths if recorded else ())
        pats = parse_patterns(patterns, reads) or (recorded.patterns if recorded else ())
        if not paths and reads:
            paths = tuple(os.path.join(project.project_dir, p) for p in DEFAULT_SEARCH_PATHS)
        if not pats and reads:
            pats = DEFAULT_PATTERNS
        if not paths or not pats:
            raise ValueError(f"Task {task_name!r}: give --search-path and --pattern key=regex")
        return Scan(tuple(os.path.abspath(p) for p in paths), pats)
    if recorded is not None:
        return recorded
    if reads:
        paths = tuple(os.path.join(project.project_dir, p) for p in DEFAULT_SEARCH_PATHS)
        return Scan(paths, DEFAULT_PATTERNS)
    return None


def scan_files(
    project: Project, task_names: Iterable[str] = (), search_paths=(), patterns=()
) -> tuple[list[str], list[str]]:
    """Write the files tables; return the messages for what was written and the problems
    (libraries without files, two libraries sharing a folder, mates missing or doubled). A
    table is written in any case, with the libraries that have files."""
    df = build_library_dataframe(project.data_set_infos, include_background=True)
    folders = {}
    for row in df.itertuples(index=False):
        folder = getattr(row, "folder_name", None)
        folders[row.library_name] = (
            folder if isinstance(folder, str) and folder else row.library_name
        )
    shared = {f: [lib for lib, g in folders.items() if g == f] for f in set(folders.values())}
    discovery = ReadDiscovery()
    written, problems = [], []
    names = scan_tasks(project, task_names)
    if (search_paths or patterns) and len(names) != 1:
        raise ValueError("--search-path and --pattern apply to one task; name it with --task")
    for name in names:
        config = project.task_configs[name]
        scan = scan_of(project, name, search_paths, patterns)
        if scan is None:
            if task_names:
                raise ValueError(
                    f"Task {name!r}: its table records no scan; give --search-path and "
                    "--pattern key=regex for its files"
                )
            continue  # other data without a recorded scan is not rescanned
        reads = config.patterns_find_reads()
        keys = list(dict.fromkeys(key for pattern in scan.patterns for key in pattern))
        rows: list[list[str]] = []
        found = 0
        for library, folder in folders.items():
            if len(shared[folder]) > 1:
                problems.append(
                    f"{name}: the libraries {', '.join(sorted(shared[folder]))} share the folder "
                    f"{folder!r}; the files of one library must be in a folder of their own"
                )
                continue
            try:
                if reads:
                    groups = discovery.find(
                        scan.search_paths, folder, scan.patterns, config.mixed_se_pe
                    )
                    if not groups:
                        raise ValueError(discovery.missing(scan.search_paths, folder))
                    for group in groups:
                        rows.append(
                            [library, group.name] + [group.paths.get(key, "") for key in keys]
                        )
                else:
                    files = discovery.find_files(scan.search_paths, folder, scan.patterns)
                    rows.append([library] + [files.get(key, "") for key in keys])
            except ValueError as e:
                problems.append(f"{name}: {library}: {e}")
                continue
            found += 1
        table = config.library_files[0]
        header = ([LIBRARY, READGROUP_COLUMN] if reads else [LIBRARY]) + keys
        _write_table(table, scan.command(name, table), header, rows)
        unit = "read groups" if reads else "rows"
        written.append(f"{name}: {len(rows)} {unit} of {found} libraries written to {table}")
    return written, problems


def _write_table(path: str, comment: str, header: list[str], rows: list[list[str]]) -> None:
    directory = os.path.dirname(os.path.abspath(path))

    def cell(value: str) -> str:
        if value and os.path.isabs(value) and value.startswith(directory + os.sep):
            return os.path.relpath(value, directory)
        return value

    os.makedirs(directory, exist_ok=True)
    with open(path, "w") as f:
        f.write(comment + "\n")
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(cell(value) for value in row) + "\n")
