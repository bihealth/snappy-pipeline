"""Abstract wrapper classes as utilities for snappy specific wrappers."""
# Note that this file tries to target a baseline of python 3.8, so outdated wrappers don't crash

import contextlib
import os
import shutil
import stat
import tempfile
import textwrap
from abc import ABCMeta, abstractmethod
from typing import Optional

from snakemake.shell import shell
from snakemake.utils import format as snakemake_format

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


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

        # Compute md5 except when filename ends with .md5
        compute_md5() {{
            fn=$1
            f=$(basename $fn)
            d=$(dirname $fn)
            pushd $d 1> /dev/null 2>&1
            md5sum $f > $f.md5
            popd 1> /dev/null 2>&1
        }}

        # Write out information about conda installation.
        conda list >{snakemake.log.conda_list}
        conda info >{snakemake.log.conda_info}
        compute_md5 {snakemake.log.conda_list}
        compute_md5 {snakemake.log.conda_info}

        # Create temp directory
        TMPDIR=$(mktemp -d)

        set -x

        # --------------------------------- Start command -----------------------------------------
    """

    footer = r"""
        # --------------------------------- End command -------------------------------------------

        set +x

        for fn in {snakemake.output}
        do
            if [[ -f "$fn" ]] && ! [[ $fn =~ \.md5$ ]]
            then
                compute_md5 $fn
            fi
        done
    """

    md5_log = r"""
        f=$(basename {log})
        d=$(dirname {log})
        pushd $d 1> /dev/null 2>&1
        md5sum $f > $f.md5
        popd 1> /dev/null 2>&1
    """

    output_links = r"""
        for path in {snakemake.output.output_links}; do
          dst=$path
          src=${{dst/\/output\//\/work\/}}
          mkdir -p "$(dirname "$dst")"
          ln -snrf "$src" "$dst"
        done
    """

    def __init__(self, snakemake, with_output_links: bool = True) -> None:
        self._snakemake = snakemake
        self._with_output_links = with_output_links
        self._check_snakemake_attributes()

    def _check_snakemake_attributes(self) -> None:
        if not getattr(self._snakemake, "log", None):
            raise AttributeError("snakemake.log is not defined")
        if not getattr(self._snakemake.log, "log", None):
            raise AttributeError("snakemake.log.log is not defined")
        if not getattr(self._snakemake.log, "conda_list", None):
            raise AttributeError("snakemake.log.conda_list is not defined")
        if not getattr(self._snakemake.log, "conda_info", None):
            raise AttributeError("snakemake.log.conda_info is not defined")

    def _create_output_links(self) -> None:
        r"""Create output/ symlinks pointing into work/ for all entries in output_links.

        Replaces the first ``/output/`` path component with ``/work/`` to
        locate the real file produced in the work directory, then creates a
        relative symlink at the output path.  Pure-Python implementation avoids
        the double-format issue that arises when the bash pattern
        ``${dst/\/output\//\/work\/}`` is processed first by Python's
        ``.format()`` and then again by Snakemake's ``shell()``.
        """
        for dst in self._snakemake.output.output_links:
            src = dst.replace("/output/", "/work/", 1)
            dst_dir = os.path.dirname(dst)
            if dst_dir:
                os.makedirs(dst_dir, exist_ok=True)
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.relpath(src, dst_dir or "."), dst)

    @abstractmethod
    def run(self, cmd: str) -> None:
        pass

    def _run(self, cmd: str, filename: Optional[str]) -> None:
        """
        Creates a temp file for the script, executes it & computes the md5 sum of the log

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
                                snakemake_format(cmd, stepout=4, snakemake=self._snakemake),
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

        shell(SnappyWrapper.md5_log.format(log=str(self._snakemake.log.log)))

        if (
            self._with_output_links
            and getattr(self._snakemake.output, "output_links", None) is not None
        ):
            self._create_output_links()


class ShellWrapper(SnappyWrapper):
    def _run_bash(self, cmd: str) -> None:
        script_log = getattr(self._snakemake.log, "script", None)
        self._run(cmd, script_log)
        if script_log:
            shell(SnappyWrapper.md5_log.format(log=script_log))

    def run(self, cmd: str) -> None:
        self._run_bash(cmd)


class RWrapper(SnappyWrapper):
    def _check_snakemake_attributes(self) -> None:
        super()._check_snakemake_attributes()
        if not getattr(self._snakemake.log, "script", None):
            raise AttributeError("snakemake.log.script is not defined")

    def _run_R(self, cmd: str) -> None:
        with open(self._snakemake.log.script, "wt") as f:
            print(snakemake_format(cmd, stepout=4, snakemake=self._snakemake), file=f)
        shell(SnappyWrapper.md5_log.format(log=self._snakemake.log.script))
        self._run(f"Rscript --vanilla {self._snakemake.log.script}", None)

    def run(self, cmd: str) -> None:
        self._run_R(cmd)


class PythonWrapper:
    """Base class for pure-Python wrappers (no external script executed).

    Re-routes ``stdout`` and ``stderr`` into the Snakemake log file declared
    via ``snakemake.log.log`` using ``contextlib.redirect_stdout`` and
    ``contextlib.redirect_stderr``.  This ensures that output produced by
    Python code -- including output from nested in-process Snakemake runs --
    ends up in the Snakemake log file regardless of whether a TTY is present.
    """

    def __init__(self, snakemake) -> None:
        self.snakemake = snakemake

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
        """Context manager capturing stdout/stderr into the Snakemake log file."""
        log_path = self._log_path()
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "at") as log_file:
            with contextlib.redirect_stdout(log_file):
                with contextlib.redirect_stderr(log_file):
                    yield

    def compute_log_md5(self) -> None:
        shell(SnappyWrapper.md5_log.format(log=self._log_path()))

    def write_conda_info(self) -> None:
        """Write conda_list/conda_info files plus their md5 sums, if declared."""
        log = getattr(self.snakemake, "log", None)
        if log is None:
            return
        for key, command in (("conda_list", "conda list"), ("conda_info", "conda info")):
            path = getattr(log, key, None)
            if not path:
                continue
            shell("{} > {}".format(command, path))
            shell(SnappyWrapper.md5_log.format(log=path))
