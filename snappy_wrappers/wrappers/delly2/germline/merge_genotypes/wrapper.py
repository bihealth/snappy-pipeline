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

    # ---------------
    # Merge genotypes
    # ---------------
    # If a single sample, there is no need to merge.
    # ``$i`` is reused from previous BCFs for-loop.
    if [[ $i -eq 1 ]]; then
        bcftools view \
            -O z \
            -o {snakemake.output.vcf} \
            $TMPDIR/cwd/1.bcf
    else
        out=$(realpath {snakemake.output.vcf})
        pushd $TMPDIR/cwd
        bcftools merge \
            -m id \
            -O z \
            -o $out \
            *.bcf
        popd
    fi
    tabix -f {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
    md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi_md5})
    """
    )
