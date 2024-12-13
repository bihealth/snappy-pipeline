"""Abstract wrapper for cnvkit.py"""

import os
import shutil
import stat
import tempfile
import textwrap

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class SimpleWrapper:
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

    def __init__(self, snakemake) -> None:
        self.snakemake = snakemake

    def run_bash(self, cmd: str) -> None:
        self._run(cmd, self.snakemake.log.sh)
        shell(SimpleWrapper.md5_log.format(log=self.snakemake.log.sh))

    def run_R(self, cmd: str) -> None:
        with open(self.snakemake.log.R, "wt") as f:
            print(cmd, file=f)
        shell(SimpleWrapper.md5_log.format(log=self.snakemake.log.R))
        self._run(f"R --vanilla < {self.snakemake.log.R}", None)

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
                            SimpleWrapper.header.format(snakemake=self.snakemake),
                            cmd,
                            SimpleWrapper.footer.format(snakemake=self.snakemake),
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

        shell(SimpleWrapper.md5_log.format(log=str(self.snakemake.log.log)))
