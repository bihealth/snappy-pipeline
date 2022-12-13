# -*- coding: utf-8 -*-
"""Wrapper for running Minimap2."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

seq_platform = snakemake.params.args["extra_infos"]["seqPlatform"]
library_kit = snakemake.params.args["extra_infos"]["libraryKit"]

shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
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
        -t 16 \
        -x $preset \
        -a {snakemake.config[step_config][ngs_mapping][minimap2][path_index]} \
        -Y \
        --MD \
        /dev/stdin \
    | samtools addreplacerg \
        -r "@RG\tID:{snakemake.wildcards.library_name}.$i\tSM:{snakemake.wildcards.library_name}\tPL:PACBIO" - \
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

# Compute MD5 sums
pushd $(dirname $out_bam)
md5sum $(basename $out_bam) >$(basename $out_bam).md5
md5sum $(basename $out_bam).bai >$(basename $out_bam).bai.md5
popd

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

# call plot-bamstats
mkdir $TMPDIR/bamstats.d
plot-bamstats \
    -p $TMPDIR/bamstats.d/ \
    {snakemake.output.report_bamstats_txt} \
|| true  # ignore failure

# Convert HTML report into one file.
inline-html \
    --in-file $TMPDIR/bamstats.d/index.html \
    --out-file {snakemake.output.report_bamstats_html} \
|| touch {snakemake.output.report_bamstats_html}

# Build MD5 files for the reports
md5sum {snakemake.output.report_bamstats_html} > {snakemake.output.report_bamstats_html_md5}
md5sum {snakemake.output.report_bamstats_txt} > {snakemake.output.report_bamstats_txt_md5}
md5sum {snakemake.output.report_flagstats_txt} >{snakemake.output.report_flagstats_txt_md5}
md5sum {snakemake.output.report_idxstats_txt} > {snakemake.output.report_idxstats_txt_md5}

# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper_bwa.py

# Logging: Save a permanent copy of the environment file used
cp $(dirname {this_file})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper_bwa.yaml
"""
)
