# -*- coding: utf-8 -*-
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

paths_tsv = " ".join(snakemake.input.tsv)
ShellWrapper(snakemake).run(
    r"""
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

# Snakemake removed the results of an earlier run (directory output)
mkdir -p {snakemake.output.calls}

gatk GermlineCNVCaller \
    --run-mode CASE \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-calls {snakemake.input.ploidy}/ploidy-calls \
    --model {snakemake.input.model}/cnv_calls-model \
    --output {snakemake.output.calls} \
    --output-prefix cnv_calls \
    {extra_arguments}
"""
)
