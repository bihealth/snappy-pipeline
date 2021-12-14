# -*- coding: utf-8 -*-
from snakemake.shell import shell

paths_tsv = " ".join(snakemake.input.tsv)
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
    --run-mode CASE \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-calls $(dirname {snakemake.input.ploidy})/ploidy-calls \
    --model {snakemake.params.args[model]} \
    --output $(dirname {snakemake.output.done}) \
    --output-prefix cnv_calls
"""
)
