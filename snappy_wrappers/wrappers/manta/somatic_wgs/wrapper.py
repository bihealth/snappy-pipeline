# -*- coding: utf-8 -*-
"""Wrapper for running Manta in somatic variant calling mode on WGS data"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

basedir=$(dirname $(dirname {snakemake.output.vcf}))
workdir=$basedir/work
outdir=$basedir/out

# Ensure the working directory is removed, configManta.py will bail out if it already exists
trap "rm -rf \"$workdir\"" EXIT

configManta.py \
    --referenceFasta {snakemake.input.reference} \
    --runDir $workdir \
    --normalBam {snakemake.input.normal_bam} \
    --tumorBam {snakemake.input.tumor_bam}

python $workdir/runWorkflow.py \
    -m local \
    -j 16

cp -ra $workdir/results $outdir
rm -rf $workdir

pushd $outdir
ln -sr results/variants/somaticSV.vcf.gz $(basename {snakemake.output.vcf})
ln -sr results/variants/somaticSV.vcf.gz.tbi $(basename {snakemake.output.vcf_tbi})
ln -sr results/variants/candidateSV.vcf.gz \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz
ln -sr results/variants/candidateSV.vcf.gz.tbi \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz.tbi

for f in results.tar.gz *.vcf.gz *.tbi; do
    md5sum $f >$f.md5
done
"""
)
