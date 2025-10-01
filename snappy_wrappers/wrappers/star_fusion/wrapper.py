# -*- coding: utf-8 -*-
"""Wrapper for running STAR-Fusion"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

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

mkdir -p output

STAR-Fusion \
    --genome_lib_dir {args[path_ctat_resource_lib]} \
    --left_fq input/reads_1.fastq.gz \
    --right_fq input/reads_2.fastq.gz \
    --output_dir output

popd
"""
)
