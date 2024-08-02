# -*- coding: utf-8 -*-
"""Wrapper for running BWA-MEM for PacBio data"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

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

i=0
for fname in $(find $(dirname {snakemake.input}) -name '*.bam' -or -name '*.fast?.gz'); do
    let "i=$i+1"

    minialign \
        -t {snakemake.config[step_config][ngs_mapping][minialign][mapping_threads]} \
        {snakemake.config[static_data_config][reference][path]} \
        $fname \
    | samtools addreplacerg \
        -r "@RT\tID:{snakemake.wildcards.library_name}.$i\tSM:{snakemake.wildcards.library_name}\tPL:PACBIO" - \
    | samtools sort -o $TMPDIR/tmp.$i.bam
done

out_bam={snakemake.output[0]}

samtools merge $out_bam $TMPDIR/tmp.*.bam
samtools index $out_bam

# Build MD5 files
pushd $(dirname {snakemake.output.bam})
md5sum $(basename {snakemake.output.bam}) > $(basename {snakemake.output.bam}).md5
md5sum $(basename {snakemake.output.bam_bai}) > $(basename {snakemake.output.bam_bai}).md5
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

# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper_bwa.py

# Logging: Save a permanent copy of the environment file used
cp $(dirname {this_file})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper_bwa.yaml
"""
)
