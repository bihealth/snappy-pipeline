# -*- coding: utf-8 -*-
"""Wrapper for computing expression values (RPKM) from counts"""

import os.path

from snakemake import shell

sample_tpl = os.path.basename(os.path.dirname(os.path.dirname(snakemake.output.tsv)))

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

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

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
cat <<"EOF" > $TMPDIR/compute_zscores.R
library( magrittr )
library( methods ) # https://github.com/tidyverse/broom/issues/67

gene_exp = read.delim("{snakemake.input.tsv}", stringsAsFactors=FALSE)

# https://stackoverflow.com/a/36159165
ranks = t(apply(-gene_exp, 1, rank))

zscores = t(apply(gene_exp, 1, scale))

sample = make.names("{sample_tpl}")
idx = which(sample == names(gene_exp))

res_df = cbind(
    rownames(gene_exp),
    zscores[ , idx],
    log2(gene_exp[ , idx] + 0.001) - log2(apply(gene_exp, 1, median) + 0.001),
    ranks[ , idx],
    gene_exp[ , idx],
    apply(gene_exp, 1, median),
    apply(gene_exp, 1, max),
    apply(gene_exp, 1, min))

colnames(res_df) = c("gene_symbol",
                     "z_score",
                     "logFC_to_cohort_median",
                     "cohort_rank",
                     "normalized_expression_in_sample",
                     "cohort_median",
                     "cohort_max",
                     "cohort_min")

write.table(res_df, file = "{snakemake.output.tsv}", quote=F, row.names=F, sep="\t")

warnings()
sessionInfo()
EOF

Rscript --vanilla $TMPDIR/compute_zscores.R
"""
)
