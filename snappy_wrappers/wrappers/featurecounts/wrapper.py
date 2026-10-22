# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

strand={args[strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(jq -r '.decision' {snakemake.input.decision})
fi

# only use primary alignments to prevent featurecounts from re-sorting the bam on disk
# re-sort BAM by query name ("unsorted" in STAR lingo) for featurecounts

bam=$(realpath {snakemake.input.bam})
pushd $TMPDIR
samtools view -h -F {args[exclude_flags]} -q {args[min_mapq]} $bam \
    | samtools sort -n -@ 2 \
    | featureCounts \
    -T {snakemake.threads} \
    -g {args[gene_attribute]} \
    -t {args[feature_type]} \
    -a {snakemake.input.features} \
    -s ${{strand}} \
    -p \
    --verbose \
    {extra_arguments} \
    -o feature_counts.tsv
popd
cp $TMPDIR/feature_counts.tsv {snakemake.output.tsv}
cp $TMPDIR/feature_counts.tsv.summary {snakemake.output.summary}
"""
)
