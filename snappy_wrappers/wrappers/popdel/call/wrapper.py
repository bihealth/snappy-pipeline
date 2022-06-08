# -*- coding: utf-8 -*-
"""Wrapper for running "popdel call"."""

import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.input.profile), file=tmpf)
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

    for name in $(cat {tmpf.name}); do
        echo $name >>$TMPDIR/profiles.txt
    done

    popdel call \
        -r {snakemake.wildcards.chrom}:{snakemake.wildcards.begin}-{snakemake.wildcards.end} \
        -o $TMPDIR/tmp.vcf \
        $TMPDIR/profiles.txt

    cat >$TMPDIR/header.txt <<EOF
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
EOF

    for sample in $(bcftools view --header-only $TMPDIR/tmp.vcf | grep '^#CHROM' | cut -f 10-); do
        echo -e "$sample\t$(echo $sample | rev | cut -d . -f 1 | rev)" >>$TMPDIR/samples.txt
    done

    bcftools annotate \
        --header $TMPDIR/header.txt \
        $TMPDIR/tmp.vcf \
    | awk -F $'\t' \
        'BEGIN {{ OFS = FS }}
         /^#/ {{ print $0; }}
         /^[^#]/ {{ $8 = $8 ";SVMETHOD=POPDELv1.1.0"; print; }}' \
    | bcftools reheader \
        --samples $TMPDIR/samples.txt \
    | bcftools sort \
        -T $TMPDIR \
        -O z \
        -o {snakemake.output.vcf}

    tabix -f {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
    """
    )
