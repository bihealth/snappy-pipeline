# -*- coding: utf-8 -*-
"""Wrapper for running PB Honey Spots on Germline data
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export LC_ALL=C

module purge
module load Python/2.7.9-foss-2015a
module load parallel/20160322-foss-2015a

PBSUITE=/fast/users/mholtgr/scratch/build/PBSuite_15.8.24/

. $PBSUITE/setup.sh

set -euo pipefail

inputs=$(echo {snakemake.input} | tr ' ' '\n' | grep '\.bam$')

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

honey_py_spots()
{{
    set -ex
    mkdir -p $1 && pushd $1

    ln -s ../ins.bam .
    ln -s ../ins.bam.bai .
    ln -s ../del.bam .
    ln -s ../del.bam.bai .

    pwd
    ls -lh
    tree

    # Call insertions
    Honey.py spots \
        --chrom $1 \
        -q 20 \
        -m 30 \
        --reference {snakemake.config[static_data_config][reference][path]} \
        -i 20 \
        -e 1 \
        -E 1 \
        --spanMax 100000 \
        --consensus None \
        ins.bam

    # Call deletions
    Honey.py spots \
        --chrom $1 \
        -q 20 \
        -m 40 \
        --reference {snakemake.config[static_data_config][reference][path]} \
        -i 20 \
        -e 2 \
        -E 2 \
        --spanMax 100000 \
        --consensus None \
        del.bam

    cat <(grep DEL del.hon.spots) <(grep INS ins.hon.spots) \
    | sort -k1,1 -k2,2g \
    >hon.spots

    cat <(grep DEL hon.spots \
          | tr ';=' '\t' \
          | awk '(($4=="DEL") && ($13 > 1) && ($5>=30) && ($13/$21)>=0.25) {{ print $0; }}' \
          | awk '{{ w=int($5); s=$3-$2; n=w<s?w:s; x=w>s?w:s; if (n/x > 0.1) {{ print $0; }} }}' \
          | awk '{{ lw=int($5/2); rw=$5-lw; mp=int(($2+$3)/2);
                   print $1 "\t" mp-lw "\t" mp+rw "\tDEL:" lw+rw ":" int($13) "/" int($21); }}') \
        <(grep INS hon.spots \
          | tr ';=' '\t' \
          | awk '(($5>=30) && ($13 > 1) && ($13/$21)>=0.25) {{ print $0; }}' \
          | awk '($4=="INS") {{
                print $1 "\t" $3 "\t" $3+1 "\t" "INS:" $5 ":" int($13) "/" int($21); }}') \
    | sort -k1,1 -k2,2g \
    >hon.bed
}}
export -f honey_py_spots

outdir=$PWD/$(dirname {snakemake.output.bed})

for input in $inputs; do
    mkdir -p $TMPDIR/$(basename $input)
    ln -sr $input $TMPDIR/$(basename $input)/ins.bam
    ln -sr $input $TMPDIR/$(basename $input)/del.bam
    ln -sr $input.bai $TMPDIR/$(basename $input)/ins.bam.bai
    ln -sr $input.bai $TMPDIR/$(basename $input)/del.bam.bai
    pushd $TMPDIR/$(basename $input)

    parallel -t -j {snakemake.config[step_config][wgs_sv_calling][pb_honey_spots][num_threads]} honey_py_spots ::: {{1..22}} X Y

    prefix=$(basename $input .bam | rev | cut -d . -f 2- | rev)
    suffix=$(basename $input .bam | rev | cut -d . -f 1 | rev)

    cat $TMPDIR/$(basename $input)/*/hon.bed \
    | sort -k1,1 -k2,2g \
    >$outdir/$prefix.pb_honey_spots.$suffix.bed

    mkdir -p $outdir/full_result
    cp -R $TMPDIR/$(basename $input) $outdir/full_result

    popd
done

pushd $(dirname {snakemake.output.bed})
for f in *.bed; do
    md5sum $f >$f.md5
done
"""
)
