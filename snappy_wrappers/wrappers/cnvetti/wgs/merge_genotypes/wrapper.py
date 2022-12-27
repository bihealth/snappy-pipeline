# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS merge_genotypes step."""

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

# Merge genotypes across cohort -------------------------------------------------------------------

bcftools merge \
    -m id \
    -O z \
    -o {snakemake.output.vcf} \
    {snakemake.input.vcf}

tabix -f {snakemake.output.vcf}

# Compute MD5 checksums ---------------------------------------------------------------------------

pushd $(dirname "{snakemake.output.vcf}")
md5sum $(basename "{snakemake.output.vcf}") >$(basename "{snakemake.output.vcf}").md5
md5sum $(basename "{snakemake.output.vcf_tbi}") >$(basename "{snakemake.output.vcf_tbi}").md5
"""
)
