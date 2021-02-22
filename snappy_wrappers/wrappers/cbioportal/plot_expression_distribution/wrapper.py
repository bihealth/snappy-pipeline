# -*- coding: utf-8 -*-
"""Wrapper for computing expression signatures from counts"""

from snakemake import shell
import os.path

if len(snakemake.output.pdf) == 1:
    snakemake.output.pdf = snakemake.output.pdf[0]
else:
    raise

sample_tpl = os.path.basename(os.path.dirname(os.path.dirname(snakemake.output.pdf)))
sample_dir = os.path.dirname(os.path.dirname(snakemake.output.pdf))

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

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
cat <<"EOF" > $TMPDIR/compute_signatures.R
library( magrittr )
library( methods ) # https://github.com/tidyverse/broom/issues/67

gene_exp = read.delim("{snakemake.input.tsv}", stringsAsFactors=FALSE)

# genes of interest
# TODO make configurable
goi = read.table(
    "/fast/groups/cubi/projects/2016-09-27_DKTK/pipeline/CLEMENS_misc/Vogelstein/Vogelstein_gene_names.txt",
    quote="\"", comment.char="")

sample = make.names("{sample_tpl}")
idx = which(sample == names(gene_exp))

m = as.matrix(gene_exp)
m = m[rownames(m) %in% goi$V1, ]

print(dim(m))

dir.create("{sample_dir}")
pdf("{snakemake.output.pdf}")

for (i in 1:nrow(m)) {{
  if (is.na(sum(m[i, ]))) continue
  plot(density(m[i, ]), main=rownames(m)[i])
  abline(v=m[i, sample], col=2, lwd=3)
}}
dev.off()

warnings()
sessionInfo()
EOF

Rscript --vanilla $TMPDIR/compute_signatures.R
"""
)
