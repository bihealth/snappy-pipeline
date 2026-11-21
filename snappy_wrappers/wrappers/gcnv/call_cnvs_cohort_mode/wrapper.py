# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
paths_tsv = " ".join(snakemake.input.tsv)

ShellWrapper(snakemake).run(
    r"""
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

mkdir -p {snakemake.output.calls}

gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L {snakemake.input.interval_list_shard} \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-calls {snakemake.input.ploidy}/ploidy-calls \
    --annotated-intervals {snakemake.input.intervals} \
    --interval-merging-rule {args[interval_merging_rule]} \
    --output {snakemake.output.calls} \
    --output-prefix cnv_calls \
    {extra_arguments}
"""
)
