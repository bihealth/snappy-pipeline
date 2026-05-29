# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
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
strand={args[strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(cat {snakemake.input.decision})
fi

# ----------------------------------------------------------------------------
# dupradar
# ----------------------------------------------------------------------------

mkdir $TMPDIR/dupradar

# Write helper script and call R
#
cat << __EOF > $TMPDIR/dupradar/run_dupradar.R
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

snake_log=${PWD}/$(dirname {snakemake.log.log})

pushd $TMPDIR/dupradar
Rscript --vanilla run_dupradar.R \
    "{snakemake.input.bam}" \
    results.tsv \
    {snakemake.input.dupradar_path_annotation_gtf} \
    ${strand} \
    ${paired_cmd} \
    {args[num_threads]} \
    "."
popd

mv $TMPDIR/dupradar/results.tsv {snakemake.output.dupradar}
"""
)
