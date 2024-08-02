# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FastQC: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

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

outdir=$(dirname $(echo {snakemake.output}  | tr ' ' '\n' | tail -n 1))

_JAVA_OPTIONS="-Xms256m -Xmx512m -XX:CompressedClassSpaceSize=512m" \
fastqc \
    --noextract \
    -o $outdir \
    -t {snakemake.params[args][num_threads]} \
    $(echo {snakemake.input} {snakemake.params[args][more_reads]} | tr ' ' '\n' | grep 'fastq.gz$\|fastq$\|sam$\|bam$')

pushd $outdir
pwd
ls -lh
for path in $(echo {snakemake.output} | tr ' ' '\n' | grep -v '.md5$' | tail -n +2); do
    fname=$(basename $path)
    md5sum $fname > $path.md5
done
"""
)
