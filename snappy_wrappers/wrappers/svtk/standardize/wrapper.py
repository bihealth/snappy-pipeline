# -*- coding: utf-8 -*-
"""Wrapper for running ``svtk normalize``"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

shell(
    r"""
    caller=$(echo {snakemake.wildcards.caller} | sed -e 's/delly2/delly/g')

    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    bcftools view \
        -O z \
        -o $TMPDIR/calls.vcf.gz \
        {snakemake.input.calls}

    svtk standardize \
        $TMPDIR/calls.vcf.gz \
        $TMPDIR/tmp.vcf \
        $caller

    bcftools sort -O z -o {snakemake.output.vcf} $TMPDIR/tmp.vcf

    tabix -s 1 -b 2 -e 2 -f {snakemake.output.vcf}
    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
