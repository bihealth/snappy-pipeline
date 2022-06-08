# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper for filtering WGS SV results for overlap with interesting regions.
"""


from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# shell.executable('/bin/bash') # XXX

if snakemake.wildcards.regions == "whole_genome":
    # We will not use path_bed below but need a value anyway.
    path_bed = "/dev/null"
else:
    path_bed = snakemake.config["step_config"]["wgs_sv_filtration"]["region_beds"][
        snakemake.wildcards.regions
    ]

shell(
    r"""
set -x

case "{snakemake.wildcards.regions}" in
    whole_genome)
        cp {snakemake.input.vcf} {snakemake.output.vcf}
        cp {snakemake.input.vcf}.tbi {snakemake.output.vcf}.tbi
    ;;
    *)
        bedtools intersect -header -wa -a {snakemake.input.vcf} -b {path_bed} \
        | sort -k1,1V -k2,2n \
        | uniq \
        | bgzip -c \
        > {snakemake.output.vcf}
        tabix -f {snakemake.output.vcf}
    ;;
esac

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
