# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's "pre-filter" and "post-filter" steps
"""

from snakemake.shell import shell

__author__ = "Nina Thiessen"
__email__ = "nina.thiessen@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo '{snakemake.params.membership}' \
| tr ' ' '\n' \
>$TMPDIR/samples.tsv

delly filter \
    --filter somatic \
    --samples $TMPDIR/samples.tsv \
    --outfile {snakemake.output.bcf} \
    {snakemake.input.bcf}

tabix -f {snakemake.output.bcf}

pushd $(dirname {snakemake.output.bcf})
md5sum $(basename {snakemake.output.bcf}) >$(basename {snakemake.output.bcf}).md5
md5sum $(basename {snakemake.output.bcf}).csi >$(basename {snakemake.output.bcf}).csi.md5
"""
)
