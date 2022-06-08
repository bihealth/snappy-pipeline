# -*- coding: utf-8 -*-
"""Wrapper for computing expression values (RPKM) from counts"""

import os.path

from snakemake import shell

rcode = os.path.join(os.path.dirname(__file__), "scripts.R")

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
cat <<"EOF" > $TMPDIR/compute_rpkm.R
library( magrittr )
library( methods ) # https://github.com/tidyverse/broom/issues/67
source("{rcode}")

mapping_df <- read.table("{snakemake.input.tsv}", sep="\t", header=1, stringsAsFactors=FALSE)
rpkm <- compute_rpkm(mapping_df)
write.table(rpkm, file="{snakemake.output.tsv}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

warnings()
sessionInfo()
EOF

Rscript --vanilla $TMPDIR/compute_rpkm.R
"""
)
