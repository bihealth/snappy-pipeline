# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py"""

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
# Gather statistics from STAR output log
# ----------------------------------------------------------------------------

get_value() {{
    fn=$1
    pattern=$2
    grep -q "$pattern" $fn
    if [ $? -ne 0 ]
    then
        echo "0"
    else
        grep "$pattern" $fn | sed -e "s/.*	//"
    fi
}}

# Heuristics to find the STAR output from the log directory
star_log=$(readlink -f {snakemake.input.bam})
star_log=$(dirname ${{star_log}})
star_log=${{star_log}}/../log/Log.final.out

if [ -r ${{star_log}} ]
then
    awk -v v=0 '{{if(v==1){{print $0}}; if(match($0, "^ *$")>0){{v=1}}}}' \
        ${{star_log}} \
        | grep -v ':$' | sed -e "s/^ *//" | sed -e "s/ |//" \
        > ${{TMPDIR}}/read_alignment_report.tsv

    n_input=$(get_value ${{TMPDIR}}/read_alignment_report.tsv "Number of input reads")
    n_unique=$(get_value ${{TMPDIR}}/read_alignment_report.tsv "Uniquely mapped reads number")
    n_multi=$(get_value ${{TMPDIR}}/read_alignment_report.tsv "Number of reads mapped to multiple loci")
    n_toomany=$(get_value ${{TMPDIR}}/read_alignment_report.tsv "Number of reads mapped to too many loci")
    n_chimera=$(get_value ${{TMPDIR}}/read_alignment_report.tsv "Number of chimeric reads")

    n_unmapped=$(echo ${{n_input}} - ${{n_unique}} - ${{n_multi}} - ${{n_toomany}} - ${{n_chimera}} | bc -lq)
    pc_unmapped=$(echo "scale=2;100*${{n_unmapped}}/${{n_input}}" | bc -ql)

    # Workaround to avoid obtaining read counts from fastq files
    n_orig=${{n_input}}

    echo "Number of unmapped reads	${{n_unmapped}}" >> ${{TMPDIR}}/read_alignment_report.tsv
    echo "Unmapped reads in %	${{pc_unmapped}}" >> ${{TMPDIR}}/read_alignment_report.tsv
    sed -i "1i\Number of original reads	${{n_orig}}" ${{TMPDIR}}/read_alignment_report.tsv
else
    touch ${{TMPDIR}}/read_alignment_report.tsv
fi

mv ${{TMPDIR}}/read_alignment_report.tsv {snakemake.output.stats}
md5sum {snakemake.output.stats} > {snakemake.output.stats_md5}
"""
)
