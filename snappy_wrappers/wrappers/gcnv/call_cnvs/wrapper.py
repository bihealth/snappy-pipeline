# -*- coding: utf-8 -*-

from snakemake.shell import shell

print("snakemake.input =", vars(snakemake.input))

shell(
    r"""
set -x
set -euo pipefail

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"

gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L {snakemake.input.interval_list_shard} \
   $(while read tsv; do echo -I $tsv; done < {snakemake.input.chunk}) \
    --contig-ploidy-calls $(dirname {snakemake.input.ploidy})/ploidy-calls \
    --annotated-intervals {snakemake.input.intervals} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output $(dirname {snakemake.output.done}) \
    --output-prefix cnv_calls
"""
)
