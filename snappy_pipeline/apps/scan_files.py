"""``snappy scan-files``: write the files tables of the ``external_data`` tasks with a ``scan``.

For every library of the sample sheets, the files are searched in a folder named like the library
(the sheet's ``folder_name``, else the library name) below the ``scan`` search paths, with the
``scan`` search patterns, as :class:`~snappy_pipeline.reads.ReadDiscovery` does. The rows go to
the first table of ``library_files``; paths below the table's directory are written relative to
it, others absolute.
"""

from __future__ import annotations

import os
from collections.abc import Iterable

from snappy_pipeline.library_table import build_library_dataframe
from snappy_pipeline.orchestration import Project
from snappy_pipeline.reads import LIBRARY, READGROUP_COLUMN, ReadDiscovery


def scan_tasks(project: Project, task_names: Iterable[str] = ()) -> list[str]:
    """Return the names of the ``external_data`` tasks with a ``scan`` section, or of those of
    ``task_names``; raises ``ValueError`` for a name without one."""
    names = [
        name
        for name, config in project.task_configs.items()
        if getattr(config, "scan", None) is not None
    ]
    if task_names:
        for name in task_names:
            if name not in names:
                raise ValueError(f"Task {name!r} is no external_data task with a scan section")
        return list(task_names)
    return names


def scan_files(project: Project, task_names: Iterable[str] = ()) -> tuple[list[str], list[str]]:
    """Write the files tables of the tasks with a ``scan``; return the messages for what was
    written and the problems (libraries without files, two libraries sharing a folder, mates
    missing or doubled). The table is written in any case, with the libraries that have files."""
    df = build_library_dataframe(project.data_set_infos, include_background=True)
    folders = {}
    for row in df.itertuples(index=False):
        folder = getattr(row, "folder_name", None)
        folders[row.library_name] = (
            folder if isinstance(folder, str) and folder else row.library_name
        )
    discovery = ReadDiscovery()
    written, problems = [], []
    for name in scan_tasks(project, task_names):
        config = project.task_configs[name]
        reads = config.patterns_find_reads()
        roots = [os.path.abspath(path) for path in config.scan.search_paths]
        patterns = config.scan.search_patterns
        keys = list(dict.fromkeys(key for pattern in patterns for key in pattern))
        shared = {f: [lib for lib, g in folders.items() if g == f] for f in set(folders.values())}
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
                    groups = discovery.find(roots, folder, patterns, config.mixed_se_pe)
                    if not groups:
                        raise ValueError(discovery.missing(roots, folder))
                    for group in groups:
                        rows.append(
                            [library, group.name] + [group.paths.get(key, "") for key in keys]
                        )
                else:
                    files = discovery.find_files(roots, folder, patterns)
                    rows.append([library] + [files.get(key, "") for key in keys])
            except ValueError as e:
                problems.append(f"{name}: {library}: {e}")
                continue
            found += 1
        table = config.library_files[0]
        _write_table(table, [LIBRARY, READGROUP_COLUMN] + keys if reads else [LIBRARY] + keys, rows)
        unit = "read groups" if reads else "rows"
        written.append(f"{name}: {len(rows)} {unit} of {found} libraries written to {table}")
    return written, problems


def _write_table(path: str, header: list[str], rows: list[list[str]]) -> None:
    directory = os.path.dirname(os.path.abspath(path))

    def cell(value: str) -> str:
        if value and os.path.isabs(value) and value.startswith(directory + os.sep):
            return os.path.relpath(value, directory)
        return value

    os.makedirs(directory, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(cell(value) for value in row) + "\n")
