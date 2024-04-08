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

strand={snakemake.config[step_config][gene_expression_quantification][strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(cat {snakemake.input.decision})
fi

# only use primary alignments to prevent featurecounts from re-sorting the bam on disk
# re-sort BAM by query name ("unsorted" in STAR lingo) for featurecounts

bam=$(realpath {snakemake.input.bam})
pushd $TMPDIR
samtools view -h -F 260 -q 255 $bam \
    | samtools sort -n -@ 2 \
    | featureCounts \
    -T 2 \
    -g gene_id \
    -t exon \
    -a {snakemake.config[step_config][gene_expression_quantification][featurecounts][path_annotation_gtf]} \
    -s ${{strand}} \
    -p \
    --verbose \
    -o feature_counts.tsv
popd
cp $TMPDIR/feature_counts.tsv {snakemake.output.tsv}
cp $TMPDIR/feature_counts.tsv.summary {snakemake.output.summary}
md5sum {snakemake.output.tsv} > {snakemake.output.tsv_md5}
md5sum {snakemake.output.summary} > {snakemake.output.summary_md5}
"""
)
