# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for LOHHLA: Snakemake wrapper.py"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Clemens Messerschmidt"

args = snakemake.params.args
normal_bam = os.path.basename(snakemake.input.normal_bam)

ShellWrapper(snakemake).run(
    r"""
# Setup subdirectories in the automatically created TMPDIR
mkdir -p $TMPDIR/out
mkdir -p $TMPDIR/bams
mkdir -p {snakemake.output.out}

# need to link bams into the same directory
ln -sr {snakemake.input.normal_bam} -t $TMPDIR/bams
ln -sr {snakemake.input.normal_bai} -t $TMPDIR/bams
ln -sr {snakemake.input.tumor_bam} -t $TMPDIR/bams
ln -sr {snakemake.input.tumor_bai} -t $TMPDIR/bams

lohhla_script=$(which lohhla)

Rscript $lohhla_script \
    --BAMDir $TMPDIR/bams/  \
    --HLAexonLoc {snakemake.input.hla_dat} \
    --HLAfastaLoc {snakemake.input.hla_fasta} \
    --cleanUp FALSE \
    --fishingStep FALSE \
    --gatkDir {snakemake.input.picard_dir} \
    --hlaPath $(realpath {snakemake.input.hla} ) \
    --mappingStep TRUE \
    --minCoverageFilter {args[min_coverage_filter]} \
    --normalBAMfile $TMPDIR/bams/{normal_bam} \
    --novoDir $(dirname $(realpath $lohhla_script)) \
    --outputDir $(realpath {snakemake.output.out}) \
    --patientId {args[tumor_library]} \
    {extra_arguments}
"""
)
