# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"

PRIORS=$TMPDIR/ploidy_priors.tsv
echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" \
> $PRIORS
for i in {{1..22}}; do
    echo -e "$i\t0\t0.01\t0.98\t0.01" >> $PRIORS
done
echo -e "X\t0.01\t0.49\t0.49\t0.01" >> $PRIORS
echo -e "Y\t0.495\t0.495\t0.01\t0" >> $PRIORS

set -x

gatk DetermineGermlineContigPloidy \
    -L {snakemake.input.interval_list} \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(while read tsv; do echo -I $tsv; done < {snakemake.input.chunk}) \
    --contig-ploidy-priors $PRIORS \
    --output $(dirname {snakemake.output}) \
    --output-prefix ploidy
"""
)
