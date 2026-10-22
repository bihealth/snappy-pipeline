# -*- coding: utf-8 -*-
"""Wrapper for running "popdel profile"."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"


args = snakemake.params.args
profile_regions = "\n".join(args["profile_regions"])

ShellWrapper(snakemake).run(
    r"""
cat >$TMPDIR/intervals.txt <<"EOF"
{profile_regions}
EOF

# The intervals have "chr" prefixes; drop them if the first contig of the reference has none
if ! head -n 1 {snakemake.input.reference_fai} | grep -q '^chr'; then
    perl -p -i -e 's/chr//g' $TMPDIR/intervals.txt
fi

OMP_NUM_THREADS=2 \
popdel profile \
    --merge-read-groups \
    -i $TMPDIR/intervals.txt \
    -o {snakemake.output.profile} \
    {snakemake.input.bam}
"""
)
