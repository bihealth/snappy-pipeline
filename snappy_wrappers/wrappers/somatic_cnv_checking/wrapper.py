# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

rscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "cnv-check-plot.R")

reference = snakemake.config["static_data_config"]["reference"]["path"]

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
# trap "rm -rf $TMPDIR" EXIT

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

md5() {{
    filename=$1
    fn=$(basename $filename)
    pushd $(dirname $filename) 1> /dev/null 2>&1
    rslt=$(md5sum $fn)
    popd 1> /dev/null 2>&1
    echo "$rslt"
}}

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5 {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5 {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
R --vanilla --slave << __EOF
source("{rscript}")

genome_lengths <- chromosome_lengths("{reference}") |> dplyr::mutate(n=dplyr::row_number())

x <- vcf_to_table("{snakemake.input.vcf}", sample="{snakemake.wildcards[library_name]}")
x <- x |> dplyr::left_join(genome_lengths, by="CHROM") |> dplyr::mutate(x=POS + Offset)

y <- read.table("{snakemake.input.tsv}", sep="\t", header=1, check.names=FALSE)
y <- y |> dplyr::mutate(Call=cn_to_call(CN))
y <- y |> dplyr::left_join(genome_lengths, by="CHROM") |> dplyr::mutate(from=start + Offset, to=stop + Offset)

pdf("{snakemake.output.cnv}", height=6.22, width=9.33)
plot_cnv(x, scale="log2") + ggplot2::ggtitle("{snakemake.wildcards[library_name]}")
plot_cnv(x, scale="sqrt") + ggplot2::ggtitle("{snakemake.wildcards[library_name]}")
dev.off()

pdf("{snakemake.output.locus}", height=6.22, width=12.44)
plot_locus(x, genome_lengths |> dplyr::filter(CHROM %in% x[["CHROM"]])) +
    ggplot2::ggtitle("{snakemake.wildcards[library_name]}")
dev.off()

pdf("{snakemake.output.segment}", height=6.22, width=12.44)
plot_segment(y, genome_lengths |> dplyr::filter(CHROM %in% y[["CHROM"]]))
dev.off()
__EOF

md5 {snakemake.output.cnv} > {snakemake.output.cnv_md5}
md5 {snakemake.output.locus} > {snakemake.output.locus_md5}
md5 {snakemake.output.segment} > {snakemake.output.segment_md5}
"""
)

shell(
    r"""
md5() {{
    filename=$1
    fn=$(basename $filename)
    pushd $(dirname $filename) 1> /dev/null 2>&1
    rslt=$(md5sum $fn)
    popd 1> /dev/null 2>&1
    echo "$rslt"
}}
md5 {snakemake.log.log} > {snakemake.log.log_md5}
"""
)
