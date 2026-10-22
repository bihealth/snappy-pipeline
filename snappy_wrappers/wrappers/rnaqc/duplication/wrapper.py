# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for RSeQC read duplication: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
# Duplication (from RSeQC package)
mkdir ${{TMPDIR}}/duplication
read_duplication.py -i {snakemake.input.bam} -o ${{TMPDIR}}/duplication/out {extra_arguments}

mv ${{TMPDIR}}/duplication/out.seq.DupRate.xls {snakemake.output.dupl_seq}
mv ${{TMPDIR}}/duplication/out.pos.DupRate.xls {snakemake.output.dupl_pos}
"""
)
