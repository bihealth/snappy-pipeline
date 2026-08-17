# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

rscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "cnv-check-plot.R")

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# Write helper script and call R
R --vanilla --slave << __EOF
source("{rscript}")

genome_lengths <- chromosome_lengths("{snakemake.input.reference}") |> dplyr::mutate(n=dplyr::row_number())

x <- vcf_to_table("{snakemake.input.vcf}", sample="{args[library_name]}")
x <- x |> dplyr::left_join(genome_lengths, by="CHROM") |> dplyr::mutate(x=POS + Offset)

y <- read.table("{snakemake.input.tsv}", sep="\t", header=1, check.names=FALSE)
y <- y |> dplyr::mutate(Call=cn_to_call(CN))
y <- y |> dplyr::mutate(CHROM=as.character(CHROM)) |> dplyr::left_join(genome_lengths, by="CHROM") |> dplyr::mutate(from=start + Offset, to=stop + Offset)

pdf("{snakemake.output.cnv}", height=6.22, width=9.33)
plot_cnv(x, scale="log2") + ggplot2::ggtitle("{args[library_name]}")
plot_cnv(x, scale="sqrt") + ggplot2::ggtitle("{args[library_name]}")
dev.off()

pdf("{snakemake.output.locus}", height=6.22, width=12.44)
plot_locus(x, genome_lengths |> dplyr::filter(CHROM %in% x[["CHROM"]])) +
    ggplot2::ggtitle("{args[library_name]}")
dev.off()

pdf("{snakemake.output.segment}", height=6.22, width=12.44)
plot_segment(y, genome_lengths |> dplyr::filter(CHROM %in% y[["CHROM"]]))
dev.off()
__EOF
"""
)
