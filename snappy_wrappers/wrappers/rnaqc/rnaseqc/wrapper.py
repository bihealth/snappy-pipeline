# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for RNA-SeQC: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
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

# RNA-SeQC analysis (coverage from 3' UTR, ...)
mkdir ${{TMPDIR}}/rnaseqc

if [[ ${{paired}} -eq 0 ]]; then
    paired_cmd=" -singleEnd "
else
    paired_cmd=""
fi

# IMPORTANT NOTE-
# The GTF annotation file (snakemake.input.rnaseqc_path_annotation_gtf)
# assumes that:
# - all records have a "transcript_id" entry among their attributes. Many records won't have it,
#   for example all "gene" (in feature column) are missing it, and it will trigger an error whn present.
# - the rRNA information is taken from the "transcript_type" attribute, following some GENCODE version.
#   For ENSEMBL, this information is in column 2 (source).

jar_path=${{JAVA_HOME}}/share/rna-seqc-1.1.8-2/RNA-SeQC_v1.1.8.jar

${{JAVA_HOME}}/bin/java -jar ${{jar_path}} \
    -r {snakemake.input.reference} \
    -t {snakemake.input.rnaseqc_path_annotation_gtf} \
    -s "Sample,{snakemake.input.bam}, " \
    ${{paired_cmd}} \
    -o ${{TMPDIR}}/rnaseqc

mv ${{TMPDIR}}/rnaseqc/metrics.tsv {snakemake.output.rnaseqc_metrics}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_low.txt {snakemake.output.rnaseqc_meancov_low}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_medium.txt {snakemake.output.rnaseqc_meancov_medium}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_high.txt {snakemake.output.rnaseqc_meancov_high}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_low.txt {snakemake.output.rnaseqc_meannorm_low}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_medium.txt {snakemake.output.rnaseqc_meannorm_medium}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_high.txt {snakemake.output.rnaseqc_meannorm_high}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_low.txt {snakemake.output.rnaseqc_gaplen_low}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_medium.txt {snakemake.output.rnaseqc_gaplen_medium}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_high.txt {snakemake.output.rnaseqc_gaplen_high}
"""
)
