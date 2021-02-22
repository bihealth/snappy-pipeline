# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for deconstructSigs
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt"

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

cat <<"EOF" >$TMPDIR/deconstructSigs.R
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
EOF

Rscript --vanilla $TMPDIR/deconstructSigs.R

md5sum {snakemake.output.tsv} > {snakemake.output.tsv}.md5
md5sum {snakemake.output.pdf} > {snakemake.output.pdf}.md5
"""
)
