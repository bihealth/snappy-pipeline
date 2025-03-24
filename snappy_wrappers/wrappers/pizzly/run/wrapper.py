# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Kallisto+Pizzly: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

shell.executable("/bin/bash")

shell(
    r"""
set -euo pipefail
set -x
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id

workdir=$(dirname {snakemake.output.done})
inputdir=$workdir/input

mkdir -p $inputdir

if [[ ! -f "$inputdir/reads_1.fastq.gz" ]]; then
    cat {args[left]} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    cat {args[right]} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

test -f output/fusion.txt \
|| kallisto quant \
    -i {args[kallisto_index]} \
    --fusion \
    -o output \
    input/reads_1.fastq.gz \
    input/reads_2.fastq.gz

mkdir -p output

pizzly \
    -k {args[kmer_size]} \
    --cache index.cache.txt \
    --align-score 2 \
    --insert-size 400 \
    --fasta {args[transcripts_fasta]} \
    --output fusions \
    --gtf {args[annotations_gtf]} \
    output/fusion.txt

pizzly_flatten_json.py fusions.json > fusions.txt

popd
"""
)
