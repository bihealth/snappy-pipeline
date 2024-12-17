"""Abstract wrapper classes as utilities for snappy specific wrappers."""

import os
import shutil
import stat
import tempfile
import textwrap
from abc import abstractmethod, ABCMeta

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class SnappyWrapper(metaclass=ABCMeta):
    header = r"""
        # Pipe everything to log file
        if [[ -n "{snakemake.log.log}" ]]; then
            if [[ "$(set +e; tty; set -e)" != "" ]]; then
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                exec &> >(tee -a "{snakemake.log.log}" >&2)
            else
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                echo "No tty, logging disabled" >"{snakemake.log.log}"
            fi
        fi

        # Compute md5
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
            if ! [[ $fn =~ \.md5$ ]]
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
          src=work/${{dst#output/}}
          ln -sr $src $dst
        done
    """

    def __init__(self, snakemake, with_output_links: bool = True) -> None:
        self._snakemake = snakemake
        self._with_output_links = with_output_links
        self._check_snakemake_attributes()

    def _check_snakemake_attributes(self):
        if not getattr(self._snakemake, "log", None):
            raise AttributeError("snakemake.log is not defined")
        if not getattr(self._snakemake.log, "log", None):
            raise AttributeError("snakemake.log.log is not defined")
        if not getattr(self._snakemake.log, "conda_list", None):
            raise AttributeError("snakemake.log.conda_list is not defined")
        if not getattr(self._snakemake.log, "conda_info", None):
            raise AttributeError("snakemake.log.conda_info is not defined")
        if not getattr(self._snakemake.log, "script", None):
            raise AttributeError("snakemake.log.script is not defined")

    @abstractmethod
    def run(self, cmd: str) -> None:
        pass

    def _run(self, cmd: str, filename: str | None) -> None:
        """
        Creates a temp file for the script, executes it & computes the md5 sum of the log

        The shell script is first created as a temporary file, and then copied over to
        the log directory.
        This allows R scripts to be saved in the log directory, rather than the uninformative
        shell script starting R.

        : param cmd: The command string (after snakemake input/output/params expansion)
        : param filename: the path where to save the script
        """
        with tempfile.NamedTemporaryFile(mode="wt", delete_on_close=False) as f:
            tempfilename = f.name

            print(
                textwrap.dedent(
                    "\n".join(
                        (
                            SnappyWrapper.header.format(snakemake=self._snakemake),
                            cmd,
                            SnappyWrapper.footer.format(snakemake=self._snakemake),
                        )
                    )
                ),
                file=f,
            )

            f.flush()
            f.close()

            current_permissions = stat.S_IMODE(os.lstat(tempfilename).st_mode)
            os.chmod(tempfilename, current_permissions | stat.S_IXUSR)

            if filename is not None:
                shutil.copy(tempfilename, filename)

            shell(tempfilename)

        shell(SnappyWrapper.md5_log.format(log=str(self._snakemake.log.log)))

        if (
            self._with_output_links
            and getattr(self._snakemake.output, "output_links", None) is not None
        ):
            shell(SnappyWrapper.output_links.format(snakemake=self._snakemake))


class ShellWrapper(SnappyWrapper):
    def _run_bash(self, cmd: str) -> None:
        self._run(cmd, self._snakemake.log.script)
        shell(SnappyWrapper.md5_log.format(log=self._snakemake.log.script))

    def run(self, cmd: str) -> None:
        self._run_bash(cmd)


class RWrapper(SnappyWrapper):
    def _run_R(self, cmd: str) -> None:
        with open(self._snakemake.log.script, "wt") as f:
            print(cmd, file=f)
        shell(SnappyWrapper.md5_log.format(log=self._snakemake.log.script))
        self._run(f"Rscript --vanilla {self._snakemake.log.script}", None)

    def run(self, cmd: str) -> None:
        self._run_R(cmd)
