# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

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
# RNA-SeQC analysis (coverage from 3' UTR, ...)
# ----------------------------------------------------------------------------

mkdir ${{TMPDIR}}/rnaseqc

if [[ ${{paired}} -eq 0 ]]; then
    paired_cmd=" -singleEnd "
else
    paired_cmd=""
fi

# IMPORTANT NOTE-
# The GTF annotation file (snakemake.config[step_config][gene_expression_quantification][rnaseqc][rnaseqc_path_annotation_gtf])
# assumes that:
# - all records have a "transcript_id" entry among their attributes. Many records won't have it,
#   for example all "gene" (in feature column) are missing it, and it will trigger an error whn present.
# - the rRNA information is taken from the "transcript_type" attribute, following some GENCODE version.
#   For ENSEMBL, this information is in column 2 (source).

jar_path=${{JAVA_HOME}}/share/rna-seqc-1.1.8-2/RNA-SeQC_v1.1.8.jar

${{JAVA_HOME}}/bin/java -jar ${{jar_path}} \
    -r {snakemake.config[static_data_config][reference][path]} \
    -t {snakemake.config[step_config][gene_expression_quantification][rnaseqc][rnaseqc_path_annotation_gtf]} \
    -s "Sample,{snakemake.input.bam}, " \
    ${{paired_cmd}} \
    -o ${{TMPDIR}}/rnaseqc

mv ${{TMPDIR}}/rnaseqc/metrics.tsv {snakemake.output.rnaseqc_metrics}
md5sum {snakemake.output.rnaseqc_metrics} > {snakemake.output.rnaseqc_metrics_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_low.txt {snakemake.output.rnaseqc_meancov_low}
md5sum {snakemake.output.rnaseqc_meancov_low} > {snakemake.output.rnaseqc_meancov_low_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_medium.txt {snakemake.output.rnaseqc_meancov_medium}
md5sum {snakemake.output.rnaseqc_meancov_medium} > {snakemake.output.rnaseqc_meancov_medium_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverage_high.txt {snakemake.output.rnaseqc_meancov_high}
md5sum {snakemake.output.rnaseqc_meancov_high} > {snakemake.output.rnaseqc_meancov_high_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_low.txt {snakemake.output.rnaseqc_meannorm_low}
md5sum {snakemake.output.rnaseqc_meannorm_low} > {snakemake.output.rnaseqc_meannorm_low_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_medium.txt {snakemake.output.rnaseqc_meannorm_medium}
md5sum {snakemake.output.rnaseqc_meannorm_medium} > {snakemake.output.rnaseqc_meannorm_medium_md5}
mv ${{TMPDIR}}/rnaseqc/meanCoverageNorm_high.txt {snakemake.output.rnaseqc_meannorm_high}
md5sum {snakemake.output.rnaseqc_meannorm_high} > {snakemake.output.rnaseqc_meannorm_high_md5}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_low.txt {snakemake.output.rnaseqc_gaplen_low}
md5sum {snakemake.output.rnaseqc_gaplen_low} > {snakemake.output.rnaseqc_gaplen_low_md5}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_medium.txt {snakemake.output.rnaseqc_gaplen_medium}
md5sum {snakemake.output.rnaseqc_gaplen_medium} > {snakemake.output.rnaseqc_gaplen_medium_md5}
mv ${{TMPDIR}}/rnaseqc/gapLengthHist_high.txt {snakemake.output.rnaseqc_gaplen_high}
md5sum {snakemake.output.rnaseqc_gaplen_high} > {snakemake.output.rnaseqc_gaplen_high_md5}
"""
)
