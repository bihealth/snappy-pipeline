# -*- coding: utf-8 -*-
"""Wrapper to create graphs of MHC alleles for HLA-LA"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

mkdir $TMPDIR

wget -O $TMPDIR/graphs.tar.gz {snakemake.params.url}
actual_md5=$(md5sum $TMPDIR/graphs.tar.gz | sed -e "s/ .*//")
if [[ $actual_md5 != {snakemake.params.expected_md5} ]]
then
    echo "Wrong checksum of {snakemake.params.url}. Expected {snakemake.params.expected_md5}, found $actual_md5"
    exit 1
fi

tar -zxvf $TMPDIR/graphs.tar.gz -C $(dirname {snakemake.output.done})

HLA-LA.pl --prepareGraph=1 \
    --customGraphDir $(dirname {snakemake.output.done})

touch {snakemake.output.done}
"""
)
