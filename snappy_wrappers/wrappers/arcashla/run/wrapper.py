# -*- coding: utf-8 -*-
"""Wrapper for running ARCAS-HLA"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

genotype_args = (
    f"--population {args['population']} --min_count {args['min_count']} "
    f"--tolerance {args['tolerance']} --max_iterations {args['max_iterations']} "
    f"--drop_threshold {args['drop_threshold']} "
    f"--zygosity_threshold {args['zygosity_threshold']}"
)
if "drop_iterations" in args:
    genotype_args += f" --drop_iterations {args['drop_iterations']}"
single_end_args = ""
if "avg" in args:
    single_end_args += f" --avg {args['avg']}"
if "std" in args:
    single_end_args += f" --std {args['std']}"

ShellWrapper(snakemake).run(
    r"""
# Single-end reads, if no read is paired in sequencing
if [[ $(samtools view -c -f 1 -@ {snakemake.threads} {snakemake.input.bam}) -eq 0 ]]; then
    single="--single"
    single_end_args="{single_end_args}"
else
    single=""
    single_end_args=""
fi

mkdir $TMPDIR/tmp
arcasHLA extract --verbose --threads {snakemake.threads} $single \
    --temp $TMPDIR/tmp \
    --outdir $TMPDIR/extracted \
    {snakemake.input.bam}

rm -rf $TMPDIR/tmp

mkdir $TMPDIR/tmp
arcasHLA genotype --verbose --threads {snakemake.threads} \
    --outdir $TMPDIR/genotyped \
    --temp $TMPDIR/tmp \
    {genotype_args} $single_end_args \
    {extra_arguments} \
    $TMPDIR/extracted/*.fq.gz

cp $TMPDIR/genotyped/*.genotype.json {snakemake.output.json}

# arcasHLA writes no allele list
touch {snakemake.output.txt}
"""
)
