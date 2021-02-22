# -*- coding: utf-8 -*-

import os
from snakemake.shell import shell

paths_tsv = " ".join(snakemake.input.tsv)

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
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-calls $(dirname {snakemake.input.ploidy})/ploidy-calls \
    --annotated-intervals {snakemake.input.intervals} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output $(dirname {snakemake.output.done}) \
    --output-prefix cnv_calls
"""
)
