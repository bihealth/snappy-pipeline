# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS coverage step."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Get preset and individual settings from configuration.
cnvetti_config = snakemake.config["step_config"]["wgs_cnv_calling"]["cnvetti"]
preset_name = cnvetti_config["preset"]
preset = cnvetti_config["presets"][preset_name]

window_length = cnvetti_config.get("window_length") or preset["window_length"]
count_kind = cnvetti_config.get("count_kind") or preset["count_kind"]
normalization = cnvetti_config.get("normalization") or preset["normalization"]

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

# Compute coverage and normalize ------------------------------------------------------------------

REF={snakemake.config[static_data_config][reference][path]}

cnvetti cmd coverage \
    -vvv \
    --considered-regions GenomeWide \
    --count-kind {count_kind} \
    --window-length {window_length} \
    --reference $REF \
    --output $TMPDIR/cov.bcf \
    --input {snakemake.input.bam}

cnvetti cmd normalize \
    -vvv \
    --normalization {normalization} \
    --input $TMPDIR/cov.bcf \
    --output {snakemake.output.bcf}

# Compute MD5 checksums ---------------------------------------------------------------------------

pushd $(dirname "{snakemake.output.bcf}")
md5sum $(basename "{snakemake.output.bcf}") >$(basename "{snakemake.output.bcf}").md5
md5sum $(basename "{snakemake.output.csi}") >$(basename "{snakemake.output.csi}").md5
"""
)
