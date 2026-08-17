# -*- coding: utf-8 -*-
"""Wrapper for running "popdel profile"."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
cat >$TMPDIR/intervals.txt <<"EOF"
chr1:35000000-36000000
chr2:174000000-175000000
chr3:36500000-37500000
chr4:88000000-89000000
chr5:38000000-39000000
chr6:38000000-39000000
chr7:38000000-39000000
chr8:19000000-20000000
chr9:19000000-20000000
chr10:19000000-20000000
chr11:19000000-20000000
chr12:19000000-20000000
chr13:25000000-26000000
chr14:25000000-26000000
chr15:25000000-26000000
chr16:25000000-26000000
chr17:31000000-32000000
chr18:31000000-32000000
chr19:31000000-32000000
chr20:33000000-34000000
chr21:21000000-22000000
chr22:25000000-26000000
EOF

if [[ "{args[reference]}" =~ .*hs?37.* ]]; then
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
