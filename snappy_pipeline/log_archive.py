# -*- coding: utf-8 -*-
"""Archive the logs of one ``snappy run`` (plans.md F3).

The Snakemake log of a run lists every job it started together with the job's log files, and
the slurm executor logs the slurm log file of every job it submits. The archive holds these
files plus a copy of the snkmt database, so building it reads one log instead of walking the
task directories.
"""

import contextlib
import re
import sqlite3
import sys
import tarfile
import tempfile
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
