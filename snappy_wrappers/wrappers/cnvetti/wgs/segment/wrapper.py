# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS segment step."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

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

# Compute segmentation ----------------------------------------------------------------------------

cnvetti cmd segment \
    -vvv \
    --segmentation {args[segmentation]} \
    --input {snakemake.input.bcf} \
    --output {snakemake.output.windows_bcf} \
    --output-segments {snakemake.output.segments_bcf}

# Compute MD5 checksums ---------------------------------------------------------------------------

pushd $(dirname "{snakemake.output.windows_bcf}")
md5sum $(basename "{snakemake.output.windows_bcf}") >$(basename "{snakemake.output.windows_bcf}").md5
md5sum $(basename "{snakemake.output.windows_bcf_csi}") >$(basename "{snakemake.output.windows_bcf_csi}").md5
md5sum $(basename "{snakemake.output.segments_bcf}") >$(basename "{snakemake.output.segments_bcf}").md5
md5sum $(basename "{snakemake.output.segments_bcf_csi}") >$(basename "{snakemake.output.segments_bcf_csi}").md5
"""
)
