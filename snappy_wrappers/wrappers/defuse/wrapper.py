# -*- coding: utf-8 -*-
"""Wrapper for running Defuse"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

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
    zcat {snakemake.params.args[left]} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    zcat {snakemake.params.args[right]} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

defuse_run.pl \
    -d {snakemake.params.args[path_dataset_directory]} \
    -1 input/reads_1.fastq.gz \
    -2 input/reads_2.fastq.gz \
    -o output \
    -p 8

popd
"""
)
