# -*- coding: utf-8 -*-
"""Wrapper for running Manta in germline variant calling mode on WGS data
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

basedir=$(dirname $(dirname {snakemake.output.vcf}))
workdir=$basedir/work
outdir=$basedir/out

# Ensure the working directory is removed, configManta.py will bail out if it already exists
trap "rm -rf \"$workdir\"" EXIT

configManta.py \
    --referenceFasta {snakemake.config[static_data_config][reference][path]} \
    --runDir $workdir \
    $(echo "{snakemake.input}" | tr ' ' '\n' | grep -v 'bai$' | sed 's/^/--bam /g')

perl -p -i -e 's/isEmail = .*/isEmail = False/g' $workdir/runWorkflow.py

python2 $workdir/runWorkflow.py \
    --jobs 16

cp -ra $workdir/results $outdir
rm -rf $workdir

pushd $outdir
tar czf results.tar.gz results
ln -sr results/variants/diploidSV.vcf.gz $(basename {snakemake.output.vcf})
ln -sr results/variants/diploidSV.vcf.gz.tbi $(basename {snakemake.output.tbi})
ln -sr results/variants/candidateSV.vcf.gz \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz
ln -sr results/variants/candidateSV.vcf.gz.tbi \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz.tbi

for f in results.tar.gz *.vcf.gz *.tbi; do
    md5sum $f >$f.md5
done
"""
)
