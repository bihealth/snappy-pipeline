# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Hera: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -euo pipefail
set -x
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id

export TMPDIR=$(mktemp -d)

workdir=$(dirname {snakemake.output.done})
inputdir=$workdir/input

mkdir -p $inputdir

if [[ ! -f "$inputdir/reads_1.fastq.gz" ]]; then
    cat {snakemake.params.args[left]} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    cat {snakemake.params.args[right]} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

hera quant \
    -i {snakemake.config[step_config][somatic_gene_fusion_calling][hera][path_index]} \
    -f {snakemake.config[step_config][somatic_gene_fusion_calling][hera][path_genome]} \
    -t 8 \
    -o $PWD \
    input/reads_1.fastq.gz \
    input/reads_2.fastq.gz

popd
"""
)
