"""Abstract wrapper classes as utilities for snappy specific wrappers."""
# Snakemake runs wrappers with the Python of their environment if it is at least 3.7
# (snakemake.script.MIN_PY_VERSION), so this file must stay valid Python 3.7.

import contextlib
import hashlib
import os
import shutil
import stat
import sys
import tempfile
import textwrap
from abc import ABCMeta, abstractmethod
from types import SimpleNamespace

from snakemake.shell import shell
from snakemake.utils import format as snakemake_format

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


def write_declared_md5s(snakemake) -> None:
    """Write ``<file>.md5`` for each output whose ``.md5`` is a declared output, too.

    Only result files declare an ``.md5`` (SODAR needs one per file). The format is that of
    ``md5sum``.
    """
    for md5_path in map(str, snakemake.output):
        if not md5_path.endswith(".md5"):
            continue
        path = md5_path[: -len(".md5")]
        digest = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                digest.update(chunk)
        if os.path.lexists(md5_path):
            os.remove(md5_path)  # may be a link to the md5 file of an input
        with open(md5_path, "wt") as f:
            f.write("{}  {}\n".format(digest.hexdigest(), os.path.basename(path)))


#: The suffixes of the scripts a job leaves next to its log: the shell script of the job and the
#: R or Python script a wrapper embeds (the log archive takes them along)
SCRIPT_SUFFIXES = ("sh", "R", "py")


def script_path(log: str, suffix: str = "sh") -> str:
    """A script of a job next to its log: ``<name>.<suffix>`` for ``<name>.log``."""
    stem = log[: -len(".log")] if log.endswith(".log") else log
    return "{}.{}".format(stem, suffix)


class SnappyWrapper(metaclass=ABCMeta):
    header = r"""
        #!/usr/bin/env bash
        set -euo pipefail

        # Pipe everything to the snakemake log file while keeping the output
        # on the original stdout (terminal or slurm job log).
        if [[ -n "{snakemake.log.log}" ]]; then
            rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
            exec > >(tee -a "{snakemake.log.log}") 2>&1
        fi

        # Temp directory of the job, which the wrapper base class removes afterwards; a rerun of
        # this script creates it again
        export TMPDIR="{tmpdir}"
        mkdir -p "$TMPDIR"

        # The job's threads for OpenMP and BLAS: Snakemake's shell() sets these to None when a
        # wrapper calls it
        export OMP_NUM_THREADS={snakemake.threads} GOTO_NUM_THREADS={snakemake.threads} \
            OPENBLAS_NUM_THREADS={snakemake.threads} MKL_NUM_THREADS={snakemake.threads} \
            VECLIB_MAXIMUM_THREADS={snakemake.threads} NUMEXPR_NUM_THREADS={snakemake.threads}

        set -x

        # --------------------------------- Start command -----------------------------------------
    """

    footer = r"""
        # --------------------------------- End command -------------------------------------------

        set +x
    """

    def __init__(self, snakemake) -> None:
        self._snakemake = snakemake
        self._check_snakemake_attributes()

    def _check_snakemake_attributes(self) -> None:
        if not getattr(self._snakemake, "log", None):
            raise AttributeError("snakemake.log is not defined")
        if not getattr(self._snakemake.log, "log", None):
            raise AttributeError("snakemake.log.log is not defined")

    @property
    def extra_arguments(self) -> str:
        """The ``extra_arguments`` of the rule's params as one string, ``{extra_arguments}`` in
        the command"""
        args = getattr(self._snakemake.params, "args", None)
        value = args.get("extra_arguments") if isinstance(args, dict) else None
        return " ".join(value) if isinstance(value, list) else ""

    @abstractmethod
    def run(self, cmd: str) -> None:
        pass

    def _run(self, cmd: str, scripts=None) -> None:
        """
        Writes the scripts next to the log, executes the shell script & writes the declared
        md5 files. The scripts stay, with every placeholder filled in: ``<name>.sh`` is the
        command as run, which ``bash`` reruns; ``scripts`` (suffix -> template) are the scripts
        the command embeds, ``<name>.R`` or ``<name>.py``, named in the command as
        ``{scripts.R}`` or ``{scripts.py}`` by their absolute paths.

        :param cmd: The command string (after snakemake input/output/params expansion)
        """
        log = str(self._snakemake.log.log)
        os.makedirs(os.path.dirname(os.path.abspath(log)), exist_ok=True)
        paths = {}
        for suffix, template in (scripts or {}).items():
            paths[suffix] = os.path.abspath(script_path(log, suffix))
            with open(paths[suffix], "wt") as f:
                print(
                    textwrap.dedent(
                        snakemake_format(
                            template,
                            stepout=4,
                            snakemake=self._snakemake,
                            extra_arguments=self.extra_arguments,
                        )
                    ),
                    file=f,
                )
        script = script_path(log)
        tmpdir = tempfile.mkdtemp()
        try:
            with open(script, "wt") as f:
                print(
                    "\n".join(
                        (
                            textwrap.dedent(
                                snakemake_format(
                                    SnappyWrapper.header,
                                    stepout=4,
                                    snakemake=self._snakemake,
                                    tmpdir=tmpdir,
                                )
                            ),
                            snakemake_format(
                                cmd,
                                stepout=4,
                                snakemake=self._snakemake,
                                extra_arguments=self.extra_arguments,
                                scripts=SimpleNamespace(**paths),
                            ),
                            textwrap.dedent(
                                snakemake_format(
                                    SnappyWrapper.footer, stepout=4, snakemake=self._snakemake
                                )
                            ),
                        )
                    ),
                    file=f,
                )
            os.chmod(script, stat.S_IMODE(os.lstat(script).st_mode) | stat.S_IXUSR)
            shell(script)
        finally:
            # The temp directory goes regardless of exceptions; the script stays
            shutil.rmtree(tmpdir, ignore_errors=True)

        write_declared_md5s(self._snakemake)


class ShellWrapper(SnappyWrapper):
    def _run_bash(self, cmd: str, scripts=None) -> None:
        self._run(cmd, scripts)

    def run(self, cmd: str, **scripts: str) -> None:
        """Run ``cmd``, a bash template. ``scripts`` are the scripts it embeds, by suffix:
        ``R=code`` or ``py=code``, written next to the log as ``<name>.R`` or ``<name>.py`` with
        their placeholders filled in, and named in ``cmd`` as ``{scripts.R}`` or ``{scripts.py}``.
        """
        self._run_bash(cmd, scripts)


class RWrapper(ShellWrapper):
    """A wrapper whose command is an R script: ``run(code)`` writes it next to the log as
    ``<name>.R``, with its placeholders filled in, and runs it with ``Rscript --vanilla``."""

    def run(self, cmd: str) -> None:
        self._run_bash("Rscript --vanilla {scripts.R}", {"R": cmd})


class PythonWrapper:
    """Base class for pure-Python wrappers (no external script executed).

    Follows the log contract of ``ShellWrapper``: ``run(func)`` calls ``func`` with stdout and
    stderr in ``snakemake.log.log``, then writes the declared md5 files.
    """

    def __init__(self, snakemake) -> None:
        self.snakemake = snakemake
        self._log_path()  # raises if snakemake.log.log is missing

    def run(self, func) -> None:
        """Run ``func()`` with its output in the log, then write the declared md5 files."""
        log_path = self._log_path()
        if os.path.exists(log_path):
            os.remove(log_path)  # the log holds this job only, as with ShellWrapper
        with self.logging_context():
            func()
        write_declared_md5s(self.snakemake)

    def _log_path(self) -> str:
        log = getattr(self.snakemake, "log", None)
        if not getattr(log, "log", None):
            raise AttributeError("snakemake.log.log is not defined")
        return os.path.realpath(str(log.log))

    @contextlib.contextmanager
    def logging_context(self):
        """Context manager capturing stdout/stderr into the Snakemake log file.

        Redirects the file descriptors as well, so the output of subprocesses (``shell()``,
        R) ends up in the log, too.
        """
        log_path = self._log_path()
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "at") as log_file:
            sys.stdout.flush()
            sys.stderr.flush()
            saved = os.dup(1), os.dup(2)
            os.dup2(log_file.fileno(), 1)
            os.dup2(log_file.fileno(), 2)
            try:
                with contextlib.redirect_stdout(log_file):
                    with contextlib.redirect_stderr(log_file):
                        yield
            finally:
                log_file.flush()
                os.dup2(saved[0], 1)
                os.dup2(saved[1], 2)
                os.close(saved[0])
                os.close(saved[1])
