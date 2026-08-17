# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py export"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

ShellWrapper(snakemake).run(
    r"""
set -x

export TMPDIR=$(mktemp -d)

cnvkit.py export bed {snakemake.input} -o $TMPDIR/out.bed
bgzip -c $TMPDIR/out.bed > {snakemake.output.bed}
tabix -f {snakemake.output.bed}

cnvkit.py export seg {snakemake.input} -o {snakemake.output.seg}

cnvkit.py export vcf {snakemake.input} -o $TMPDIR/out.vcf
bgzip -c $TMPDIR/out.vcf > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

rm -rf $TMPDIR
"""
)
