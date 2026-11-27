# -*- coding: utf-8 -*-
"""Archive the logs of one ``snappy run``.

The Snakemake log of a run lists every job it started together with the job's log files, and
the slurm executor logs the slurm log file of every job it submits. The archive holds these
files plus a copy of the snkmt database, so building it reads one log instead of walking the
task directories.

The start hook of a run writes the record ``<run start>.git.txt`` next to the Snakemake log: the
version of snappy-pipeline with the commit it runs from, and, in a git repository, the commit of
the configuration and its uncommitted changes. The archive holds that record too.
"""

import contextlib
import importlib.metadata
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path

import snappy_pipeline
from snappy_pipeline import __version__

#: ``log:`` line of a job in the Snakemake log, with comma-separated paths
_JOB_LOGS = re.compile(r"^    log: (.+)$", re.M)

#: Message of the slurm executor for a submitted job
_SLURM_LOG = re.compile(r"has been submitted with SLURM jobid \S+ \(log: (.+)\)\.$", re.M)

#: snkmt database below the project directory, as passed by ``snappy run``
SNKMT_DB = Path(".snakemake", "log", "snkmt.sqlite")


def run_log_files(snakemake_log: Path) -> list[str]:
    """Return the job and slurm log files that the Snakemake log of a run names."""
    text = Path(snakemake_log).read_text(errors="replace")
    files = [path.strip() for line in _JOB_LOGS.findall(text) for path in line.split(",")]
    files += _SLURM_LOG.findall(text)
    return list(dict.fromkeys(path for path in files if path))


@dataclass(frozen=True)
class GitState:
    """The commit of the project directory and what differs from it in the tracked files."""

    commit: str
    #: ``git status --porcelain`` lines of the tracked files, empty when clean
    status: str
    diff: str

    @property
    def changed_files(self) -> list[str]:
        return [line[3:] for line in self.status.splitlines()]

    @property
    def text(self) -> str:
        record = f"commit {self.commit}\n"
        if self.status:
            record += f"uncommitted changes:\n{self.status}\n{self.diff}"
        return record


def run_git(directory: Path | str, *args: str) -> subprocess.CompletedProcess | None:
    """Run git in ``directory`` and return the result, ``None`` without git."""
    git = shutil.which("git")
    if git is None:
        return None
    return subprocess.run([git, "-C", str(directory), *args], capture_output=True, text=True)


def git_state(directory: Path) -> GitState | None:
    """Return the git state of ``directory``, or ``None`` without git, outside a repository or
    before its first commit."""
    head = run_git(directory, "rev-parse", "HEAD")
    if head is None or head.returncode != 0:
        return None
    status = run_git(directory, "status", "--porcelain", "--untracked-files=no").stdout
    diff = run_git(directory, "diff", "HEAD").stdout if status else ""
    return GitState(head.stdout.strip(), status, diff)


def checkout_commit(directory: Path) -> tuple[str, bool] | None:
    """The commit of the git checkout at ``directory`` and whether it has uncommitted changes.

    None without git, or when ``directory`` is not the top level of a repository, as a package
    installed into an environment below a repository is not that repository's checkout.
    """
    top = run_git(directory, "rev-parse", "--show-toplevel")
    if top is None or top.returncode != 0 or Path(top.stdout.strip()).resolve() != directory:
        return None
    head = run_git(directory, "rev-parse", "HEAD")
    if head is None or head.returncode != 0:
        return None
    status = run_git(directory, "status", "--porcelain", "--untracked-files=no")
    return head.stdout.strip(), bool(status and status.stdout.strip())


def direct_url_commit(direct_url: str | None) -> str | None:
    """The commit in the ``direct_url.json`` of a package installed from git, else None."""
    if not direct_url:
        return None
    return json.loads(direct_url).get("vcs_info", {}).get("commit_id") or None


def pipeline_commit() -> tuple[str, bool] | None:
    """The commit snappy-pipeline runs from and whether its checkout has uncommitted changes.

    From the checkout of an editable install, else from the ``direct_url.json`` of a package
    installed from git (``pip install git+...``); None for a wheel or sdist install.
    """
    checkout = Path(snappy_pipeline.__file__).resolve().parent.parent
    if commit := checkout_commit(checkout):
        return commit
    try:
        direct_url = importlib.metadata.distribution("snappy-pipeline").read_text("direct_url.json")
    except importlib.metadata.PackageNotFoundError:
        return None
    commit_id = direct_url_commit(direct_url)
    return (commit_id, False) if commit_id else None


def pipeline_version(full_commit: bool = False) -> str:
    """``snappy-pipeline <version>`` with the commit it runs from, when known.

    ``modified`` marks a checkout with uncommitted changes; the commit is abbreviated to 12
    characters unless ``full_commit``.
    """
    commit = pipeline_commit()
    if commit is None:
        return f"snappy-pipeline {__version__}"
    sha, modified = commit
    sha = sha if full_commit else sha[:12]
    return f"snappy-pipeline {__version__} (commit {sha}{', modified' if modified else ''})"


def untracked_sheet_hints(directory: str, tables: list[str], sheets: list[str]) -> list[str]:
    """Hints for the files tables that the run record misses: ignored by git while a sample sheet
    is tracked, or allowed but not committed; nothing outside a repository."""
    inside = run_git(directory, "rev-parse", "--is-inside-work-tree")
    if inside is None or inside.returncode != 0:
        return []

    def tracked(path: str) -> bool:
        return run_git(directory, "ls-files", "--error-unmatch", "--", path).returncode == 0

    sheet_tracked = any(tracked(sheet) for sheet in sheets)
    hints = []
    for table in tables:
        if tracked(table):
            continue
        name = os.path.relpath(table, directory)
        if run_git(directory, "check-ignore", "-q", "--", table).returncode == 0:
            if sheet_tracked:
                lines = f"`!/{name}`"
                if os.sep in name:
                    lines = f"`!/{name.split(os.sep)[0]}/` and {lines}"
                hints.append(
                    f"{name} is ignored by git while the sample sheet is tracked; add {lines} to "
                    ".gitignore and commit it, so that the run record covers the files of the "
                    "libraries"
                )
        else:
            hints.append(f"{name} is not committed; commit it with the configuration")
    return hints


def git_record_path(snakemake_log: Path) -> Path:
    """Return ``<run start>.git.txt`` next to the Snakemake log."""
    snakemake_log = Path(snakemake_log)
    return snakemake_log.with_name(snakemake_log.name.removesuffix(".snakemake.log") + ".git.txt")


def record_git_state(directory: Path, snakemake_logs: list[str]) -> None:
    """Write the record of the run next to its Snakemake log and warn about uncommitted changes;
    the start hook of the workflow.

    The record names the version and commit of snappy-pipeline and, in a git repository, the
    commit of the configuration with the diff of its uncommitted changes.
    """
    state = git_state(Path(directory))
    record = pipeline_version(full_commit=True) + "\n" + (state.text if state else "")
    for snakemake_log in snakemake_logs:
        git_record_path(snakemake_log).write_text(record)
    if state is not None and state.status:
        print(
            "[snappy] Warning: uncommitted changes in "
            + ", ".join(state.changed_files)
            + f"; the run records commit {state.commit[:12]} and the diff in .snakemake/log/",
            file=sys.stderr,
            flush=True,
        )


def archive_path(directory: Path, snakemake_log: Path) -> Path:
    """Return ``logs/<run start>.tar.gz`` below ``directory`` for the given Snakemake log."""
    stamp = Path(snakemake_log).name.removesuffix(".snakemake.log")
    return Path(directory) / "logs" / f"{stamp}.tar.gz"


def build_archive(directory: Path, snakemake_log: Path) -> Path:
    """Write the log archive of the run with ``snakemake_log`` and return its path."""
    directory, snakemake_log = Path(directory).resolve(), Path(snakemake_log).resolve()
    target = archive_path(directory, snakemake_log)
    target.parent.mkdir(exist_ok=True)

    def arcname(path: Path) -> str:
        path = path.resolve()
        return str(path.relative_to(directory)) if path.is_relative_to(directory) else path.name

    from snappy_wrappers.snappy_wrapper import SCRIPT_SUFFIXES, script_path

    with tarfile.open(target, "w:gz") as tar:
        tar.add(snakemake_log, arcname=arcname(snakemake_log))
        git_record = git_record_path(snakemake_log)
        if git_record.is_file():
            tar.add(git_record, arcname=arcname(git_record))
        for name in run_log_files(snakemake_log):
            # A job log and the scripts of the job next to it
            scripts = [directory / script_path(name, suffix) for suffix in SCRIPT_SUFFIXES]
            for path in [directory / name, *scripts]:
                if path.is_file():
                    tar.add(path, arcname=arcname(path))
        db = directory / SNKMT_DB
        if db.is_file():
            # snkmt may still write, so copy with SQLite's backup API instead of the file
            with tempfile.TemporaryDirectory() as tmp:
                copy = Path(tmp) / db.name
                with contextlib.closing(sqlite3.connect(db)) as src:
                    with contextlib.closing(sqlite3.connect(copy)) as dst:
                        src.backup(dst)
                tar.add(copy, arcname=str(SNKMT_DB))
    return target


def latest_snakemake_log(directory: Path) -> Path | None:
    """Return the newest Snakemake log of the project, or ``None``."""
    logs = sorted(Path(directory, ".snakemake", "log").glob("*.snakemake.log"))
    return logs[-1] if logs else None


def archive_after_run(directory: Path, snakemake_logs: list[str]) -> None:
    """Archive the logs from a workflow hook; report failures without failing the run."""
    for snakemake_log in snakemake_logs:
        try:
            print(f"Logs of this run: {build_archive(directory, Path(snakemake_log))}")
        except (OSError, sqlite3.Error, tarfile.TarError) as e:
            print(f"Could not archive the logs of this run: {e}", file=sys.stderr)
