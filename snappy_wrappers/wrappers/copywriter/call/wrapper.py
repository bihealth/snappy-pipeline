# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

import os

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

shell.executable("/bin/bash")

rscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "snappy-copywriter-call.R")
dest = os.path.dirname(os.path.dirname(snakemake.output.bins_txt))

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
cat <<"EOF" > $TMPDIR/run_copywriter_call.R
source("{rscript}")

plot_genes <- read.table(
    "{snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][plot_genes]}",
    sep="\t", header=1, stringsAsFactors=FALSE )

postProcess(
    destination.folder="{dest}",
    donorID="{snakemake.wildcards.library_name}",
    fullID="{snakemake.wildcards.library_name}",
    cghcall.params=list(
        inter=c(-0.4, 0.4), nclass=5, divide=5, robustsig="yes", cellularity=0.6, ncpus=4
    ),
    plot.params=list(
        genes=plot_genes,
        main="{snakemake.wildcards.library_name}",
        ylim=c( -3, 4 )
    ),
    genomeRelease="{snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][genome]}",
    features={snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][features]}
)

warnings()
EOF

cat $TMPDIR/run_copywriter_call.R

Rscript --vanilla $TMPDIR/run_copywriter_call.R
"""
)

shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
md5sum {snakemake.output.bins_txt} >{snakemake.output.bins_txt_md5}
md5sum {snakemake.output.gene_call_txt} >{snakemake.output.gene_call_txt_md5}
md5sum {snakemake.output.gene_log2_txt} >{snakemake.output.gene_log2_txt_md5}
md5sum {snakemake.output.segments_txt} >{snakemake.output.segments_txt_md5}
"""
)
