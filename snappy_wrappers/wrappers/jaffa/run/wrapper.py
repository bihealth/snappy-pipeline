# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for JAFFA: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id

# Java fun1
export MALLOC_ARENA_MAX=4

export JAFFA_REF_BASE=/fast/projects/fusionbench/software/JAFFA-1.08/JAFFA_REFERENCE_FILES_HG38_GENCODE22

workdir=$(dirname {snakemake.output.done})
inputdir=$workdir/input

mkdir -p $inputdir

if [[ ! -f "$$inputdir/reads_1.fastq.gz" ]]; then
    cat {snakemake.params.args[left]} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$$inputdir/reads_2.fastq.gz" ]]; then
    cat {snakemake.params.args[right]} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir \
&& jaffa-hybrid input/reads_1.fastq.gz input/reads_2.fastq.gz \
&& popd
"""
)
