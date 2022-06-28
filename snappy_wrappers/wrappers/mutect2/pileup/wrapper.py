# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GetPileupSummaries: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]
common_biallelic = snakemake.config["step_config"]["somatic_variant_calling"]["mutect2"][
    "common_biallelic"
]

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
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

out_base=$TMPDIR/out/$(basename {snakemake.output.pileup} .pileup)

gatk --java-options '-Xms4000m -Xmx8000m' GetPileupSummaries \
    --input {snakemake.input.bam} \
    --reference {reference} \
    --variant {common_biallelic} \
    --intervals {common_biallelic} \
    --output $out_base.pileup

pushd $TMPDIR && \
    for f in $out_base.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mkdir -p $(dirname {snakemake.output.pileup})
mv $out_base.* $(dirname {snakemake.output.pileup})
"""
)
