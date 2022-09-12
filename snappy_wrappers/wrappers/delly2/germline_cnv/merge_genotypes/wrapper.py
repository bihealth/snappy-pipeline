# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's calls tep
"""

import tempfile

from snakemake.shell import shell

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
    exec &> >(tee -a "{snakemake.log.log}")
    set -x
    # -----------------------------------------------------------------------------

    # Write out information about conda installation
    conda list > {snakemake.log.conda_list}
    conda info > {snakemake.log.conda_info}
    md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
    md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}

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
        cp $TMPDIR/cwd/1.bcf {snakemake.output.bcf}
        cp $TMPDIR/cwd/1.bcf.csi {snakemake.output.csi}
    else
        out=$(realpath {snakemake.output.bcf})
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
    md5sum $(basename {snakemake.output.bcf}) > $(basename {snakemake.output.bcf_md5})
    md5sum $(basename {snakemake.output.csi}) > $(basename {snakemake.output.csi_md5})
    """
    )

    # Compute MD5 sums of logs
    shell(
        r"""
    md5sum {snakemake.log.log} > {snakemake.log.log_md5}
    """
    )
