# -*- coding: utf-8 -*-
"""Wrapper for running Minimap2."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

args = getattr(snakemake.params, "args", {})
seq_platform = args["extra_infos"]["seqPlatform"]
library_kit = args["extra_infos"]["libraryKit"]

ShellWrapper(snakemake).run(
    r"""
set -x
mkdir -p $TMPDIR/{{out,sorted,sort.tmp}}

if [[ "{library_kit}" == "PacBio HiFi" ]]; then
    preset=map-hifi
elif [[ "{library_kit}" == "PacBio CLR" ]]; then
    preset=map-pb
elif [[ "{library_kit}" == ONT* ]]; then
    preset=map-ont
else
    >&2 echo "Unknown library kit {library_kit}"
    exit 1
fi

i=1
for fname in $(find $(dirname {snakemake.input}) -name '*.bam' -or -name '*.fast?.gz'); do
    basename=$(basename $fname .bam)

    if [[ "$fname" == *.bam ]]; then \
        samtools fastq -F 2048 $fname; \
    else \
        zcat $fname; \
    fi \
    | minimap2 \
        -t {args[mapping_threads]} \
        -x $preset \
        -a {args[path_index]} \
        -Y \
        --MD \
        /dev/stdin \
    | samtools addreplacerg \
        -r "@RG\tID:{args[library_name]}.$i\tSM:{args[library_name]}\tPL:PACBIO" - \
    >$TMPDIR/out/$i.bam

    samtools sort -m 4G -@ 3 \
        -O BAM \
        -o $TMPDIR/sorted/$i.bam \
        -T $TMPDIR/sort.tmp/ \
        $TMPDIR/out/$i.bam

    let "i=$i+1"
done

out_bam={snakemake.output.bam}

if [[ $i == 2 ]]; then
    mv $TMPDIR/sorted/1.bam $out_bam
else
    samtools merge -@ 8 $out_bam $TMPDIR/sorted/*.bam
fi

cp $out_bam /tmp

samtools index $out_bam

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}
"""
)
