# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's calls tep
"""

import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.input), file=tmpf)
    tmpf.flush()
    # Actually run the script.
    shell(
        r"""
    # -----------------------------------------------------------------------------
    # Redirect stderr to log file by default and enable printing executed commands
    exec &> >(tee -a "{snakemake.log}")
    set -x
    # -----------------------------------------------------------------------------

    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    mkdir $TMPDIR/cwd

    i=0
    for x in $(cat {tmpf.name}); do
        let "i=$i+1"
        ln -s $(readlink -f $x) $TMPDIR/cwd/$i.bcf
        ln -s $(readlink -f $x).csi $TMPDIR/cwd/$i.bcf.csi
    done

    out=$(realpath {snakemake.output.bcf})
    
    # ---------------
    # Merge genotypes
    # ---------------
    # If a single sample, there is no need to merge.
    # ``$i`` is reused from previous BCFs for-loop.
    if [[ $i -eq 1 ]]; then
        cp $TMPDIR/cwd/1.bcf {snakemake.output.bcf}
        cp $TMPDIR/cwd/1.bcf.csi {snakemake.output.bcf}".csi"
    else
        pushd $TMPDIR/cwd
        bcftools merge \
            -m id \
            -O b \
            -o $out \
            *.bcf
        popd
        tabix -f {snakemake.output.bcf}
    fi

    pushd $(dirname {snakemake.output.bcf})
    md5sum $(basename {snakemake.output.bcf}) >$(basename {snakemake.output.bcf}).md5
    md5sum $(basename {snakemake.output.bcf}).csi >$(basename {snakemake.output.bcf}).csi.md5
    """
    )
