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

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

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

# Build MD5 files for the reports
md5sum {snakemake.output.report_bamstats_txt} > {snakemake.output.report_bamstats_txt_md5}
md5sum {snakemake.output.report_flagstats_txt} >{snakemake.output.report_flagstats_txt_md5}
md5sum {snakemake.output.report_idxstats_txt} > {snakemake.output.report_idxstats_txt_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
