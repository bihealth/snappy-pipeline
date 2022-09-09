# -*- coding: utf-8 -*-
"""Wrapper for computing expression z-scores for cbioportal
given CNV gene calls and expression counts"""

import os

from snakemake import shell

rcode = os.path.join(os.path.dirname(__file__), "scripts.R")

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

rscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "snappy-run-zscores.R")

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
source( "{rscript}" )

zscores_mapping_df <- read.table("{snakemake.input.tsv}", sep="\t", head=1, stringsAsFactors=FALSE)

zscores = compute_z_scores(df = zscores_mapping_df, min_nb_diploid=5)

mapping_df <- read.table("{snakemake.input.tsv}", sep="\t", header=1, stringsAsFactors=FALSE)
zscores = compute_z_scores(df = mapping_df, min_nb_diploid=5)
write.table(zscores, file="{snakemake.output.tsv}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

warnings()
sessionInfo()
EOF

Rscript --vanilla $TMPDIR/compute_zscores.R
"""
)
