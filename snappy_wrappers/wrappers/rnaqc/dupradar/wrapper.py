# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -euo pipefail
set -x

# Setup auto-cleaned TMPDIR
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

# ----------------------------------------------------------------------------
# Inititalisation: paired & strandedness decisions
# ----------------------------------------------------------------------------
# Find out single or paired ended
n_pair=$(samtools view -f 0x1 {snakemake.input.bam} | head -n 1000 | wc -l || true)
if [[ $n_pair -eq 0 ]]; then
    paired=0
else
    paired=1
fi

# Find out strand
strand={snakemake.config[step_config][gene_expression_quantification][strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(cat {snakemake.input.decision})
fi

# ----------------------------------------------------------------------------
# dupradar
# ----------------------------------------------------------------------------

mkdir ${{TMPDIR}}/dupradar

# Write helper script and call R
#
cat << __EOF > ${{TMPDIR}}/dupradar/run_dupradar.R
library(dupRadar)

args = commandArgs(trailingOnly=TRUE)

file = args[1]
out = args[2]
gtf = args[3]
stranded = as.integer(args[4])
paired = as.logical(args[5])
threads = as.integer(args[6])
outdir = args[7]

dm = analyzeDuprates(file, gtf, stranded, paired, threads)

write.table(dm, file=out, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
__EOF

if [[ ${{paired}} -eq 0 ]]; then
    paired_cmd="FALSE"
else
    paired_cmd="TRUE"
fi

snake_log=${{PWD}}/$(dirname {snakemake.log})

pushd ${{TMPDIR}}/dupradar
Rscript --vanilla run_dupradar.R \
    "{snakemake.input.bam}" \
    results.tsv \
    {snakemake.config[step_config][gene_expression_quantification][dupradar][dupradar_path_annotation_gtf]} \
    ${{strand}} \
    ${{paired_cmd}} \
    {snakemake.config[step_config][gene_expression_quantification][dupradar][num_threads]} \
    "."
popd

mv ${{TMPDIR}}/dupradar/results.tsv {snakemake.output.dupradar}
md5sum {snakemake.output.dupradar} > {snakemake.output.dupradar_md5}
"""
)
