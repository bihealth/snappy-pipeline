# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

from snakemake import shell

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
cat << __EOF > $TMPDIR/run_dupradar.R
library(dupRadar)

args = commandArgs(trailingOnly=TRUE)

file = args[1]
out = args[2]
gtf = args[3]
stranded = as.integer(args[4])
paired = as.logical(args[5])
threads = as.integer(args[6])
outdir = args[7]
sample = basename(dirname(outdir))

dm = analyzeDuprates(file, gtf, stranded, paired, threads, tmpDir="$TMPDIR", autosort=FALSE)

bitmap("{snakemake.output.densplot}", type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
duprateExpDensPlot(DupMat=dm, main=sample)
dev.off()

bitmap("{snakemake.output.boxplot}", type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
par(mar=c(10, 4, 4, 2) + 0.1)
duprateExpBoxplot(DupMat=dm, main=sample)
dev.off()

bitmap("{snakemake.output.plot}", type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
duprateExpPlot(DupMat=dm, main=sample)
dev.off()

bitmap("{snakemake.output.hist}", type='png16m', height=1000, width=1000, unit='px', taa=4, gaa=4)
expressionHist(dm)
dev.off()

write.table(dm, file=out, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
__EOF

out_dir=$(readlink -f $(dirname {snakemake.output.tsv}))

# Find if the bam file contains paired reads
n_pair=$(samtools view -f 0x1 {snakemake.input.bam} | head -n 1000 | wc -l || true)
if [[ $n_pair -eq 0 ]]; then
    paired="FALSE"
else
    paired="TRUE"
fi

Rscript --vanilla $TMPDIR/run_dupradar.R \
    "{snakemake.input.bam}" \
    "{snakemake.output.tsv}" \
    {snakemake.config[step_config][gene_expression_quantification][dupradar][path_annotation_gtf]} \
    {snakemake.config[step_config][gene_expression_quantification][dupradar][strandedness]} \
    $paired \
    {snakemake.config[step_config][gene_expression_quantification][dupradar][num_threads]} \
    $out_dir

for f in $(ls $out_dir); do
    md5sum $out_dir/$f > $out_dir/$f.md5
done
"""
)
