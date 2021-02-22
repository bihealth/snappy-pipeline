# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GATK HC GVCF CombineGVCF: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

intervals = " ".join(snakemake.params["args"]["intervals"])

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Size of batches for joining, ensures that open file handles and GATK threads
# as well as memory usage remain limited
BATCH_SIZE=50
# Break bands at multiple of this number, must be smaller than window size and
# overlap size
BREAK_BANDS_LENGTH=10000

# TODO: add through shell.prefix
export TMPDIR=$HOME/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
#trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# TODO: Fix link problems of tabix. What the link?
export JAVA_HOME=$(dirname $(which gatk_nonfree))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

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

# Split input files into batches of $BATCH_SIZE
split_input()
{{
    iteration=$1
    shift

    outdir=$TMPDIR/result.$iteration.d
    mkdir -p $outdir

    splitdir=$TMPDIR/splitting/split.$iteration.d
    mkdir -p $splitdir

    echo $* | tr ' ' '\n' | split -a 6 -l $BATCH_SIZE - $splitdir/list.

    for list in $splitdir/list.*; do
        token=$(echo $list | rev | cut -d . -f 1 | rev)
        outname=$outdir/tmp.$iteration.$token.g.vcf.gz

        1>&2 gatk_nonfree -Xmx6g -Djava.io.tmpdir=$TMPDIR \
            --analysis_type CombineGVCFs \
            -R {snakemake.config[static_data_config][reference][path]} \
            --dbsnp {snakemake.config[static_data_config][dbsnp][path]} \
            $(while read fname; do
                echo --variant $fname
              done <$list) \
            $(for itv in {intervals}; do
                echo --intervals $itv
              done) \
            --breakBandsAtMultiplesOf $BREAK_BANDS_LENGTH \
            -o $outname
    done

    echo $outdir/tmp.$iteration.*.g.vcf.gz | sort
}}

files="{snakemake.input}"
iteration=1
num_vcf=$(echo $files | wc -w)

# Merge gVCF files, in batches
while [[ $num_vcf -gt 1 ]]; do
    files=$(split_input $iteration $files)
    num_vcf=$(echo $files | wc -w)
    let "iteration=$iteration+1"
done

# Copy out resulting file
cp $files {snakemake.output.vcf}

# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
