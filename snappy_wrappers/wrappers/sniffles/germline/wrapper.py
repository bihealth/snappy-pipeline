# -*- coding: utf-8 -*-
"""Wrapper for running Sniffles on germline data
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

module purge
module load BCFtools/1.3.1-foss-2015a

set -euo pipefail

inputs=$(echo {snakemake.input} | tr ' ' '\n' | grep '\.bam$')
outfile=$(echo {snakemake.output} | tr ' ' '\n' | grep 'vcf.gz' | head -n 1)
outdir=$(dirname $outfile)

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

for input in $inputs; do
    mkdir -p $outdir/more
    outname=$outdir/more/$(basename $input .bam).vcf

    sniffles \
        -m $input \
        -t {snakemake.config[step_config][wgs_sv_calling][sniffles][num_threads]} \
        -v $outname \
        --tmp_file $TMPDIR/sniffles.tmp \
        --cluster \
        --genotype
done

touch {snakemake.output}
"""
)
