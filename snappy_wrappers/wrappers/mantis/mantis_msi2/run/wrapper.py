# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MANTIS: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Clemens Messerschmidt"

args = getattr(snakemake.params, "args", {})


ShellWrapper(snakemake).run(
    r"""
set -x

mkdir -p $TMPDIR/out

# The default thresholds are those that the authors recommend for WES,
# but should also work for reasonable deep WGS according to them.
# https://github.com/OSU-SRLab/MANTIS/issues/25
#
# Only one thread is used, see https://github.com/OSU-SRLab/MANTIS/issues/57

mantis-msi2 \
    -t {snakemake.input.tumor_bam}  \
    -n {snakemake.input.normal_bam} \
    --genome {snakemake.input.reference} \
    --bedfile {snakemake.input.loci_bed} \
    --min-read-length {args[min_read_length]} \
    --min-read-quality {args[min_read_quality]} \
    --min-locus-quality {args[min_locus_quality]} \
    --min-locus-coverage {args[min_locus_coverage]} \
    --min-repeat-reads {args[min_repeat_reads]} \
    --threads 1 \
    -o $TMPDIR/out/$(basename {snakemake.output.result}) \
    {extra_arguments}


mv $TMPDIR/out/* $(dirname {snakemake.output.result})
"""
)

