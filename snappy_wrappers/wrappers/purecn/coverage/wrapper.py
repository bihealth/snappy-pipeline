# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing PureCN coverage"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

ShellWrapper(snakemake).run(
    r"""
# Create coverage
mkdir -p $(dirname {snakemake.output.coverage})
Rscript /opt/PureCN/Coverage.R --force \
    --seed {config[seed]} \
    --out-dir $(dirname {snakemake.output.coverage}) \
    --bam {snakemake.input.bam} \
    --intervals {snakemake.input.intervals}

# Rename coverage file name
d=$(dirname {snakemake.output.coverage})
bam_name=$(basename {snakemake.input.bam} .bam)
fn="$d/${{bam_name}}_coverage_loess.txt.gz"

test -e $fn
mv $fn {snakemake.output.coverage}
"""
)
