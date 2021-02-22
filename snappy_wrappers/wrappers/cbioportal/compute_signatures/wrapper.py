# -*- coding: utf-8 -*-
"""Wrapper for computing expression signatures from counts"""

from snakemake import shell
import os.path

# import pprint
# pprint.pprint(snakemake.output.pdf)

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

sample = make.names("{sample_tpl}")
idx = which(sample == names(gene_exp))

m = as.matrix(gene_exp)
sig_list  = expression.signatures::get_signature_list(id = "Gene_symbol")

ssGSEA = GSVA::gsva(expr = m,
                    gset.idx.list = sig_list,
                    method = "gsva", no.bootstraps = 0,
                    #method = "ssgsea", ssgsea.norm=T,
                    rnaseq = T,
                    parallel.sz = 0,
                    min.sz = 5)

sigs = ssGSEA$es.obs

dir.create("{sample_dir}")
pdf("{snakemake.output.pdf}")

for (i in 1:nrow(sigs)) {{
  signature = rownames(sigs)[i]
  plot(density(sigs[i, ]), main=rownames(sigs)[i])
  abline( v= sigs[i, sample], col=2, lwd=3)
}}
dev.off()

warnings()
sessionInfo()
EOF

Rscript --vanilla $TMPDIR/compute_signatures.R
"""
)
