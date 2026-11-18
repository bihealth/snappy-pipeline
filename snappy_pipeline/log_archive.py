# -*- coding: utf-8 -*-
"""Archive the logs of one ``snappy run``.

The Snakemake log of a run lists every job it started together with the job's log files, and
the slurm executor logs the slurm log file of every job it submits. The archive holds these
files plus a copy of the snkmt database, so building it reads one log instead of walking the
task directories.

In a git repository, the start hook of a run writes the commit of the configuration and its
uncommitted changes next to the Snakemake log (``<run start>.git.txt``), and the archive holds
that record too.
"""

import contextlib
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path

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


def git_state(directory: Path) -> GitState | None:
    """Return the git state of ``directory``, or ``None`` without git, outside a repository or
    before its first commit."""
    git = shutil.which("git")
    if git is None:
        return None

    def run(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run([git, "-C", str(directory), *args], capture_output=True, text=True)

    head = run("rev-parse", "HEAD")
    if head.returncode != 0:
        return None
    status = run("status", "--porcelain", "--untracked-files=no").stdout
    diff = run("diff", "HEAD").stdout if status else ""
    return GitState(head.stdout.strip(), status, diff)


def git_record_path(snakemake_log: Path) -> Path:
    """Return ``<run start>.git.txt`` next to the Snakemake log."""
    snakemake_log = Path(snakemake_log)
    return snakemake_log.with_name(snakemake_log.name.removesuffix(".snakemake.log") + ".git.txt")


def record_git_state(directory: Path, snakemake_logs: list[str]) -> None:
    """Write the git record of the run next to its Snakemake log and warn about uncommitted
    changes; the start hook of the workflow."""
    state = git_state(Path(directory))
    if state is None:
        return
    for snakemake_log in snakemake_logs:
        git_record_path(snakemake_log).write_text(state.text)
    if state.status:
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

    with tarfile.open(target, "w:gz") as tar:
        tar.add(snakemake_log, arcname=arcname(snakemake_log))
        git_record = git_record_path(snakemake_log)
        if git_record.is_file():
            tar.add(git_record, arcname=arcname(git_record))
        for name in run_log_files(snakemake_log):
            path = directory / name
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
