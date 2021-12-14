# -*- coding: utf-8 -*-

from snakemake.shell import shell

paths_tsv = " ".join(snakemake.input.tsv)

shell(
    r"""
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"

set -x

gatk DetermineGermlineContigPloidy \
    --model {snakemake.params.args[model]} \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --output $(dirname {snakemake.output}) \
    --output-prefix ploidy
"""
)
