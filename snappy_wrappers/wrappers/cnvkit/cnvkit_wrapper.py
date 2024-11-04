"""Abstract wrapper for cnvkit.py"""

import textwrap

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class CnvkitWrapper:
    header = r"""
        # Also pipe everything to log file
        if [[ -n "{snakemake.log.log}" ]]; then
            if [[ "$(set +e; tty; set -e)" != "" ]]; then
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                exec &> >(tee -a "{snakemake.log.log}" >&2)
            else
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                echo "No tty, logging disabled" >"{snakemake.log.log}"
            fi
        fi

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

        with open(self.snakemake.log.sh, "wt") as f:
            print(
                textwrap.dedent(
                    "\n".join((CnvkitWrapper.header, self.command, CnvkitWrapper.footer))
                ),
                file=f,
            )

        shell(self.snakemake.log.sh)

        shell(CnvkitWrapper.md5_log.format(log=str(self.snakemake.log.log)))
