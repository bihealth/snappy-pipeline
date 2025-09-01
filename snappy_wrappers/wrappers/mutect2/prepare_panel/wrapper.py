# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

intervals = getattr(snakemake.input, "intervals", "")
if intervals:
    intervals = f"--intervals {intervals}"

shell(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

vcf=$(basename --suffix=.gz {snakemake.output.vcf})

gatk --java-options "-Xmx6g" Mutect2 \
    --tmp-dir ${{tmpdir}} \
    --reference {snakemake.input.reference} \
    --input {snakemake.input.normal_bam} \
    {intervals} \
    --max-mnp-distance 0 \
    --output $tmpdir/$vcf

bgzip $tmpdir/$vcf
tabix $tmpdir/$vcf.gz

pushd $tmpdir
md5sum $vcf.gz > $vcf.gz.md5
md5sum $vcf.gz.tbi > $vcf.gz.tbi.md5
popd

mv $tmpdir/$vcf.gz $tmpdir/$vcf.gz.md5 $tmpdir/$vcf.gz.tbi $tmpdir/$vcf.gz.tbi.md5 $(dirname {snakemake.output.vcf})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
