# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for filtering to exons/regions of interest."""

import textwrap

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

if snakemake.wildcards.exon_list == "genome_wide":
    shell(
        textwrap.dedent(
            r"""
    set -x

    cp {snakemake.input.vcf} {snakemake.output.vcf}
    cp {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    """
        )
    )
else:
    bed_path = snakemake.config["step_config"]["somatic_variant_filtration"]["exon_lists"][
        snakemake.wildcards.exon_list
    ]
    shell(
        textwrap.dedent(
            r"""
    set -x

    bedtools intersect \
        -a {snakemake.input.vcf} \
        -b {bed_path} \
        -wa \
        -u \
        -g {snakemake.config[static_data_config][reference][path]}.genome \
        -sorted \
        -header \
    | bgzip -c \
    > {snakemake.output.vcf}

    tabix -f {snakemake.output.vcf}
    """
        )
    )

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.output.vcf} > {snakemake.output.vcf}.md5
md5sum {snakemake.output.vcf_tbi} > {snakemake.output.vcf_tbi}.md5
md5sum {snakemake.log} >{snakemake.log}.md5
"""
)
