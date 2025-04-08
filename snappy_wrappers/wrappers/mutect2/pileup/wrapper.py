# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GetPileupSummaries: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.input.reference
common_variants = snakemake.input.common_variants

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
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
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

out_base=$tmpdir/$(basename {snakemake.output.pileup} .pileup)

gatk --java-options '-Xms4000m -Xmx8000m' GetPileupSummaries \
    --input {snakemake.input.bam} \
    --reference {reference} \
    --variant {common_variants} \
    --intervals {common_variants} \
    --output $out_base.pileup

pushd $tmpdir && \
    for f in $out_base.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv $out_base.* $(dirname {snakemake.output.pileup})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
