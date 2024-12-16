"""Abstract wrapper for cnvkit.py"""

import os
import stat
import textwrap

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class CnvkitWrapper:
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
            ext="${{f##*.}}"
            if [[ $ext != "md5" ]]
            then
                pushd $d 1> /dev/null 2>&1
                md5sum $f > $f.md5
                popd 1> /dev/null 2>&1
            fi
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

        for fn in {snakemake.output}
        do
            compute_md5 $fn
        done
        compute_md5 {snakemake.log.sh}
    """

    md5_log = r"""
        f=$(basename {log})
        d=$(dirname {log})
        pushd $d 1> /dev/null 2>&1
        md5sum $f > $f.md5
        popd 1> /dev/null 2>&1
    """

    def __init__(self, snakemake, command) -> None:
        self.snakemake = snakemake
        self.command = command

    def preamble(self):
        pass

    def run(self) -> None:
        self.preamble()

        cmd_path = self.snakemake.log.sh
        with open(cmd_path, "wt") as f:
            print(
                textwrap.dedent(
                    "\n".join(
                        (
                            CnvkitWrapper.header.format(snakemake=self.snakemake),
                            self.command,
                            CnvkitWrapper.footer.format(snakemake=self.snakemake),
                        )
                    )
                ),
                file=f,
            )
        current_permissions = stat.S_IMODE(os.lstat(cmd_path).st_mode)
        os.chmod(cmd_path, current_permissions | stat.S_IXUSR)

        shell(cmd_path)

        shell(CnvkitWrapper.md5_log.format(log=str(self.snakemake.log.log)))
