# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS tumor_normal_ratio step.

When a matched normal sample is given for the tumor then a log2-transformed ratio is computed,
otherwise the log2-transformed relative coverage of the tumor is forwarded.
"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

tumor_bcf = snakemake.input.tumor_bcf
normal_bcf = getattr(snakemake.input, "normal_bcf", None)

ShellWrapper(snakemake).run(
    r"""
set -x

if [[ "{normal_bcf}" == "None" ]]; then
    cp {tumor_bcf} {snakemake.output.bcf}
    cp {tumor_bcf}.csi {snakemake.output.bcf}.csi
else
    tumor=$(bcftools view -h {tumor_bcf} | grep '^#CHROM' | rev | cut -f 1 | rev)
    normal=$(bcftools view -h {normal_bcf} | grep '^#CHROM' | rev | cut -f 1 | rev)

    cnvetti cmd merge-cov \
        --output $TMPDIR/merged.bcf \
        {tumor_bcf} \
        {normal_bcf}

    cnvetti cmd ratio \
        --numerator-sample $tumor \
        --denominator-sample $normal \
        --output {snakemake.output.bcf} \
        $TMPDIR/merged.bcf

    tabix -f {snakemake.output.bcf}
fi
"""
)
