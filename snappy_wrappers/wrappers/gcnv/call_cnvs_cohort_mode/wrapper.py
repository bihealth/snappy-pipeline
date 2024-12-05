# -*- coding: utf-8 -*-

from snakemake.shell import shell

paths_tsv = " ".join(snakemake.input.tsv)

print("snakemake.input =", vars(snakemake.input))

shell(
    r"""
set -x
set -euo pipefail

# We use /tmp as the temporary directory as we will need to compile things and
# this is small but needs many small file accesses. Cluster files sytems are not
# so good here but local disks are fine.
export TMPDIR=/tmp
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L {snakemake.input.interval_list_shard} \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-calls $(dirname {snakemake.input.ploidy})/ploidy-calls \
    --annotated-intervals {snakemake.input.intervals} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output $(dirname {snakemake.output.done}) \
    --output-prefix cnv_calls
"""
)
