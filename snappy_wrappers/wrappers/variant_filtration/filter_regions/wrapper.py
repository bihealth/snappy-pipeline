# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for regions filter for variant_filtration."""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

# Short-circuit in case of performing no filtration
if snakemake.wildcards.regions == "whole_genome":
    shell(
        r"""
    # Regions set to "whole_genome", just link out the data.
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
    """
    )
    sys.exit(0)  # everything went well!


# Actual Filtration -------------------------------------------------------------------------------

if snakemake.wildcards.regions == "whole_genome":
    path_bed = "/dev/null"
else:
    path_bed = snakemake.config["step_config"]["variant_filtration"]["region_beds"][
        snakemake.wildcards.regions
    ]

shell(
    r"""
set -x

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

if [[ "{snakemake.wildcards.regions}" != whole_genome ]]; then
    bedtools intersect -u -header -wa -a {snakemake.input.vcf} -b {path_bed} \
    | bcftools norm --remove-duplicates \
    | bcftools sort -o {snakemake.output.vcf} -O z
else  # else, "all"
    link=1
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
fi

if [[ "${{link-0}}" -ne 1 ]]; then
    tabix -f {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
fi
"""
)
