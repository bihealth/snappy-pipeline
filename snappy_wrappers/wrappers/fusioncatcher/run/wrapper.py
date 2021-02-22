# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Fusioncatcher: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell(
    r"""
set -ex
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id

export MALLOC_ARENA_MAX=4

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

EXE=/fast/projects/fusionbench/2016-05-18_gene_fusion_benchmark/software/fusioncatcher-2017-08-02/fusioncatcher/bin/fusioncatcher
DATA={snakemake.config[step_config][somatic_gene_fusion_calling][fusioncatcher][data_dir]}

$EXE \
    --no-update-check \
    --keep \
    --data=$DATA \
    --input=$(echo {snakemake.params.args[tumor]} | tr ' ' ',') \
    $(if [[ -n "{snakemake.params.args[normal]}" ]]; then \
        echo "--normal={snakemake.params.args[normal]}" | tr ' ' ','; \
    fi) \
    --tmp=$TMPDIR \
    --output=$(dirname {snakemake.output.done}) \
    --threads={snakemake.config[step_config][somatic_gene_fusion_calling][fusioncatcher][num_threads]}
"""
)
