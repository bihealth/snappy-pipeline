# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for region filter for wgs_mei_filtration.
"""


from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

if snakemake.wildcards.regions == "whole_genome":
    path_bed = "/dev/null"
else:
    path_bed = snakemake.config["step_config"]["wgs_mei_filtration"]["region_beds"][
        snakemake.wildcards.regions
    ]

shell(
    r"""
set -x

if [[ "{snakemake.wildcards.regions}" != whole_genome ]]; then
    bedtools intersect -header -wa -a {snakemake.input.vcf} -b {path_bed} \
    | sort -k1,1V -k2,2n \
    | uniq \
    | bgzip -c \
    > {snakemake.output.vcf}
else  # else, "all"
    cp {snakemake.input.vcf} {snakemake.output.vcf}
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
