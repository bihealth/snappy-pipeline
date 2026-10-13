# -*- coding: utf-8 -*-
"""Publish the result files of a task into its ``output/`` directory.

``work/`` holds every file a task writes and is what other tasks read. ``output/`` only holds
relative symlinks to the published files, created by one local job per task. That job declares a
single output, ``output/.manifest.tsv``, instead of one output or one job per link.
"""

import os
from pathlib import Path

#: Manifest of the published files, below the task's ``output/``
MANIFEST = ".manifest.tsv"


def work_paths(result_files: list[str]) -> list[str]:
    """Return the ``work/`` files behind a task's result files (paths below ``output/``)."""
    paths = []
    for path in result_files:
        if not path.startswith("output/"):
            raise ValueError(f"result file outside output/: {path}")
        paths.append("work/" + path.removeprefix("output/"))
    return paths


def _size(path: Path) -> int:
    if path.is_dir():
        return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    return path.stat().st_size


def _published(relative: list[Path]) -> list[Path]:
    """Return what to link: the files, but for a ``.done`` marker that is the only published
    file in its directory the directory, whose files are no declared outputs (trimmed reads,
    the result directories of the fusion callers)."""
    parents = [rel.parent for rel in relative]
    return [
        rel.parent if rel.name == ".done" and sum(p == rel.parent for p in parents) == 1 else rel
        for rel in relative
    ]


def publish(task_dir: str, work_files: list[str], manifest: str) -> None:
    """Link ``work_files`` (below ``task_dir/work``) into ``task_dir/output`` and write the
    manifest. Symlinks in ``output/`` that point to no published file are removed; other files
    there are left alone.
    """
    work_dir, output_dir = Path(task_dir, "work"), Path(task_dir, "output")
    relative = _published([Path(path).relative_to(work_dir) for path in work_files])
    links = {output_dir / rel for rel in relative}
    for dirpath, dirnames, filenames in os.walk(output_dir):
        for name in dirnames + filenames:
            path = Path(dirpath, name)
            if path.is_symlink() and path not in links:
                path.unlink()
    for rel in relative:
        link = output_dir / rel
        link.parent.mkdir(parents=True, exist_ok=True)
        target = os.path.relpath(work_dir / rel, link.parent)
        if link.is_symlink():
            if os.readlink(link) == target:
                continue
            link.unlink()
        elif link.is_dir() and not any(link.iterdir()):
            link.rmdir()  # left from linking the files of a published directory one by one
        link.symlink_to(target)
    published = {str(rel) for rel in relative}
    with open(manifest, "wt") as f:
        print("path\tsize\tmd5", file=f)
        for rel in sorted(relative):
            md5 = f"{rel}.md5" if f"{rel}.md5" in published else ""
            print(f"{rel}\t{_size(work_dir / rel)}\t{md5}", file=f)
