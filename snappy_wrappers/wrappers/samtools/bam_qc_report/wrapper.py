# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Samtools - BAM QC report: Snakemake wrapper.py"""

from snakemake import shell

shell.executable("/bin/bash")


shell(
    r"""
set -x

# Write out information about conda installation.
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}

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

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# Validate input
if [[ "{snakemake.params.args[bam_count]}" -eq  0 ]]; then
    echo "No BAM files provided!"
    exit 1
elif [[ "{snakemake.params.args[bam_count]}" -gt  1 ]]; then
    echo "Multiple BAM files provided!"
    echo "{snakemake.params.args[bam]}"
    exit 1
fi

# QC Report
samtools stats    {snakemake.params.args[bam]} > {snakemake.output.bamstats}
samtools flagstat {snakemake.params.args[bam]} > {snakemake.output.flagstats}
samtools idxstats {snakemake.params.args[bam]} > {snakemake.output.idxstats}


# Build MD5 files for the reports
md5sum {snakemake.output.bamstats}  > {snakemake.output.bamstats_md5}
md5sum {snakemake.output.flagstats} > {snakemake.output.flagstats_md5}
md5sum {snakemake.output.idxstats}  > {snakemake.output.idxstats_md5}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
"""
)
