# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for RSeQC read duplication: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# Find out single or paired ended
n_pair=$(samtools view -f 0x1 {snakemake.input.bam} | head -n 1000 | wc -l || true)
if [[ $n_pair -eq 0 ]]; then
    paired=0
else
    paired=1
fi

# Find out strand
strand={args[strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(cat {snakemake.input.decision})
fi

# Duplication (from RSeQC package)
mkdir ${{TMPDIR}}/duplication
read_duplication.py -i {snakemake.input.bam} -o ${{TMPDIR}}/duplication/out

mv ${{TMPDIR}}/duplication/out.seq.DupRate.xls {snakemake.output.dupl_seq}
mv ${{TMPDIR}}/duplication/out.pos.DupRate.xls {snakemake.output.dupl_pos}
"""
)
