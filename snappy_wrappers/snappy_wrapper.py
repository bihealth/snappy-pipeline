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
from typing import Optional

from snakemake.shell import shell
from snakemake.utils import format as snakemake_format

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


def write_declared_md5s(snakemake) -> None:
    """Write ``<file>.md5`` for each output whose ``.md5`` is a declared output, too.

    Only published files declare an ``.md5`` (SODAR needs one per file). The format is that of
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

        # Create temp directory
        TMPDIR=$(mktemp -d)

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

    def _run(self, cmd: str, filename: Optional[str]) -> None:
        """
        Creates a temp file for the script, executes it & writes the declared md5 files

        The shell script is first created as a temporary file, and then copied over to
        the log directory.
        This allows R scripts to be saved in the log directory, rather than the uninformative
        shell script starting R.

        :param cmd: The command string (after snakemake input/output/params expansion)
        :param filename: the path where to save the script
        """
        tempfilename = None
        try:
            # delete=False is safe on all Python versions.
            # It ensures the file is not unlinked upon closing the context manager.
            with tempfile.NamedTemporaryFile(mode="wt", delete=False) as f:
                tempfilename = f.name

                print(
                    textwrap.dedent(
                        "\n".join(
                            (
                                snakemake_format(
                                    SnappyWrapper.header,
                                    stepout=4,
                                    snakemake=self._snakemake,
                                ),
                                snakemake_format(
                                    cmd,
                                    stepout=4,
                                    snakemake=self._snakemake,
                                    extra_arguments=self.extra_arguments,
                                ),
                                snakemake_format(
                                    SnappyWrapper.footer,
                                    stepout=4,
                                    snakemake=self._snakemake,
                                ),
                            )
                        )
                    ),
                    file=f,
                )
                f.flush()
                # Exiting the 'with' context manager safely closes the file.

            # Since the file is closed, we can reliably adjust permissions, copy, and run it.
            current_permissions = stat.S_IMODE(os.lstat(tempfilename).st_mode)
            os.chmod(tempfilename, current_permissions | stat.S_IXUSR)

            if filename is not None:
                shutil.copy(tempfilename, filename)

            shell(tempfilename)

        finally:
            # Manually clean up the file on exit, regardless of exceptions
            if tempfilename is not None:
                try:
                    os.unlink(tempfilename)
                except OSError:
                    pass

        write_declared_md5s(self._snakemake)


class ShellWrapper(SnappyWrapper):
    def _run_bash(self, cmd: str) -> None:
        self._run(cmd, getattr(self._snakemake.log, "script", None))

    def run(self, cmd: str) -> None:
        self._run_bash(cmd)


class RWrapper(SnappyWrapper):
    def _check_snakemake_attributes(self) -> None:
        super()._check_snakemake_attributes()
        if not getattr(self._snakemake.log, "script", None):
            raise AttributeError("snakemake.log.script is not defined")

    def _run_R(self, cmd: str) -> None:
        with open(self._snakemake.log.script, "wt") as f:
            script = snakemake_format(
                cmd, stepout=4, snakemake=self._snakemake, extra_arguments=self.extra_arguments
            )
            print(script, file=f)
        self._run(f"Rscript --vanilla {self._snakemake.log.script}", None)

    def run(self, cmd: str) -> None:
        self._run_R(cmd)


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
        if log is not None:
            if getattr(log, "log", None):
                return os.path.realpath(str(log.log))
            if str(log):
                return os.path.realpath(str(log))
        raise AttributeError("snakemake.log.log is not defined")

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
