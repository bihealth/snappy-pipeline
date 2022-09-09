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

# Run rseqc to infer strandedness
infer_experiment.py \
    -r "{snakemake.config[step_config][gene_expression_quantification][strandedness][path_exon_bed]}" \
    -i "{snakemake.input.bam}" \
    > "{snakemake.output.tsv}"

md5sum {snakemake.output.tsv} > {snakemake.output.tsv_md5}

# Set strandedness based on input parameter or inferred value
strand={snakemake.config[step_config][gene_expression_quantification][strand]}
if [ ${{strand}} -eq -1 ]
then
    pattern="^This is (Pair|Single)End Data$"
    endedness=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

    pattern="^Fraction of reads failed to determine: *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    failed=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

    if [ "$endedness" == "Single" ]
    then
        pattern="^Fraction of reads explained by \"\\+\\+,\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
        pattern="^Fraction of reads explained by \"\\+\\-,\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    else
        pattern="^Fraction of reads explained by \"1\\+\\+,1\\-\\-,2\\+\\-,2\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
        pattern="^Fraction of reads explained by \"1\\+\\-,1\\-\\+,2\\+\\+,2\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    fi

    reverse=$(echo $reverse | tr '[Dde]' 'E')

    strand=0
    if [ $(echo "$forward > {snakemake.config[step_config][gene_expression_quantification][strandedness][threshold]}" | bc -l) -gt 0 ]
    then
        strand=1
    fi
    if [ $(echo "$reverse > {snakemake.config[step_config][gene_expression_quantification][strandedness][threshold]}" | bc -l) -gt 0 ]
    then
        strand=2
    fi
fi

# Write strandedness
echo ${{strand}} > {snakemake.output.decision}
md5sum {snakemake.output.decision} > {snakemake.output.decision_md5}
"""
)
