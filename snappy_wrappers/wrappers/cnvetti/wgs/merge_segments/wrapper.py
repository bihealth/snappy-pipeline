# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS merge_segments step."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

# Write out information about conda installation --------------------------------------------------

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Setup auto-cleaned TMPDIR -----------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Merge segmentation results ----------------------------------------------------------------------

cnvetti cmd merge-seg \
    -vvv \
    --output {snakemake.output.bcf} \
    {snakemake.input.bcf}

# Compute MD5 checksums ---------------------------------------------------------------------------

pushd $(dirname "{snakemake.output.bcf}")
md5sum $(basename "{snakemake.output.bcf}") >$(basename "{snakemake.output.bcf}").md5
md5sum $(basename "{snakemake.output.csi}") >$(basename "{snakemake.output.csi}").md5
"""
)
