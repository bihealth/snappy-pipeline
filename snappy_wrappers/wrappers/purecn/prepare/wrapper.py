# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing exome kit intervals for PureCN"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

mappability = getattr(snakemake.input, "mappability", "")
reptiming = getattr(snakemake.input, "reptiming", "")

ShellWrapper(snakemake).run(
    r"""
# Create panel
uncompressed=$(echo "{snakemake.output.optimized}" | sed -e "s/\.gz$//")

mkdir -p $(dirname {snakemake.output.intervals})
mkdir -p $(dirname $uncompressed)

Rscript /opt/PureCN/IntervalFile.R \
    --out-file {snakemake.output.intervals} \
    --export $uncompressed \
    --in-file {snakemake.input.bait_regions} \
    --fasta {snakemake.input.reference} --genome {config[genome_name]} \
    $(if [[ -n "{mappability}" ]]; then
        echo "--mappability {mappability}"
    fi) \
    $(if [[ -n "{reptiming}" ]]; then
        echo "--reptiming {reptiming}"
    fi)

bgzip $uncompressed
tabix {snakemake.output.optimized}
"""
)
