# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# TODO: add through shell.prefix
export TMPDIR=/fast/users/$USER/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/out

out_base=$TMPDIR/out/$(basename {snakemake.output.vcf} .vcf.gz)

gatk -Xmx16g GenomicsDBImport \
    --tmp-dir ${{TMPDIR}} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --input {snakemake.input.normal_bam} \
    --output ${{out_base}}.vcf.gz \

tabix -f ${{out_base}}.vcf.gz

rm -f ${{out_base}}.vcf.idx

tabix -f ${{out_base}}.vcf.gz

pushd $TMPDIR && \
    for f in ${{out_base}}.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv ${{out_base}}.* $(dirname {snakemake.output.vcf})
"""
)
