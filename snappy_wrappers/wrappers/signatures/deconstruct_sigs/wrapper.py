# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for deconstructSigs"""

from snappy_wrappers.snappy_wrapper import RWrapper

__author__ = "Clemens Messerschmidt"

RWrapper(snakemake).run(
    r"""
library(tidyverse)
library(deconstructSigs)

mutsigtable <- read_delim("{snakemake.input.tsv}",
  "\t", escape_double = FALSE, col_names = FALSE,
  col_types = cols(X2 = col_character()),
  trim_ws = TRUE)

colnames(mutsigtable) = c("Sample", "chr", "pos", "ref", "alt")

m = as.data.frame(mutsigtable)
m$chr = as.factor(m$chr)

sigs.input = mut.to.sigs.input(mut.ref = m)
                               #bsg = BSgenome.Hsapiens.UCSC.hg38::Hsapiens)

output.sigs = whichSignatures(tumor.ref = sigs.input,
                              #signatures.ref = signatures.nature2013,
                              signatures.ref = signatures.cosmic,
                              contexts.needed = TRUE,
                              #tri.counts.method = "exome2genome")
                              tri.counts.method = "default")

sigs = data.frame(colnames(output.sigs$weights), t(output.sigs$weights))
write_tsv(sigs, path = "{snakemake.output.tsv}")

pdf("{snakemake.output.pdf}", 7, 7)
plotSignatures(output.sigs)
dev.off()
"""
)
