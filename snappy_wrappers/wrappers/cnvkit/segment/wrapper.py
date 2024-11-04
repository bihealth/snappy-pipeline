# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py segment"""

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

if "variants" in snakemake.input:
    variants = r"""
        ---vcf {snakemake.input.variants} \
        {snakemake.params.sample_id} {snakemake.params.normal_id} \
        {snakemake.params.min_variant_depth} {snakemake.params.zygocity_freq}
    """.format(
        snakemake=snakemake,
    )
else:
    variants = ""

cmd = r"""
cnvkit.py segment --processes {snakemake.params.proceses} \
    -o {snakemake.output.segments} --dataframe {snakemake.output.dataframe} \
    --method {snakemake.params.method} --threshold {snakemake.params.threshold} {smooth_cbs} \
    {drop_low_coverage} --drop-outliers {snakemake.params.drop_outliers} \
    {variants} \
    {snakemake.input.coverage}
""".format(
    snakemake=snakemake,
    smooth_cbs="--smooth-cbs" if snakemake.params.smooth_cbs else "",
    drop_low_coverage="--drop-low-coverage" if snakemake.params.drop_low_coverage else "",
    variants=variants,
)

CnvkitWrapper(snakemake, cmd).run()
