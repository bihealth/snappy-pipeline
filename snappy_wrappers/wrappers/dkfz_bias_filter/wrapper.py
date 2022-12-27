# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for oxo-G flagging.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -euo pipefail

set -x

# "Local" TMPDIR as the scripts try to do "rename" across file sests otherwise
export TMPDIR=$(dirname $(dirname {snakemake.output.vcf}))/tmp
mkdir -p $TMPDIR
trap "rm -rf $TMPDIR" EXIT KILL TERM INT HUP

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

out={snakemake.output.vcf}

dkfzbiasfilter.py \
    --tempFolder $TMPDIR \
    --writeQC \
    {snakemake.input.vcf} \
    {snakemake.input.bam} \
    {snakemake.config[static_data_config][reference][path]} \
    ${{out%.gz}}

if [[ ! -s ${{out%.gz}} ]]; then
    bcftools view --header-only {snakemake.input.vcf} \
    > ${{out%.gz}}
fi

bgzip ${{out%.gz}}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5 && \
    popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
