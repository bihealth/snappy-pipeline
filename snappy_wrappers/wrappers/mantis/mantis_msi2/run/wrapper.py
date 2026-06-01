# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MANTIS: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Clemens Messerschmidt"

args = getattr(snakemake.params, "args", {})


ShellWrapper(snakemake).run(
    r"""
set -x

mkdir -p $TMPDIR/out

# The following config is recommended by the authors for WES,
# but should also work for reasonable deep WGS according to them.
# https://github.com/OSU-SRLab/MANTIS/issues/25
#
# Only one thread is used, see https://github.com/OSU-SRLab/MANTIS/issues/57

mantis-msi2 \
    -t {snakemake.input.tumor_bam}  \
    -n {snakemake.input.normal_bam} \
    --genome {snakemake.input.reference} \
    --bedfile {snakemake.input.loci_bed} \
    --min-read-length 35 \
    --min-read-quality 20.0 \
    --min-locus-quality 25.0 \
    --min-locus-coverage 20 \
    --min-repeat-reads 1 \
    --threads 1 \
    -o $TMPDIR/out/$(basename {snakemake.output.result})


mv $TMPDIR/out/* $(dirname {snakemake.output.result})
"""
)

