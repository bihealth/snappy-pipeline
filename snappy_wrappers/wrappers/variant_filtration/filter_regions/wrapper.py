# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for regions filter for variant_filtration."""

import os
import sys
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

args = getattr(snakemake.params, "args", {})

# Short-circuit in case of performing no filtration
if args["filter_mode"] == "whole_genome":
    from snakemake.shell import shell
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

if args["filter_mode"] == "whole_genome":
    path_bed = "/dev/null"
else:
    path_bed = args["filter_config"][args["filter_mode"]]

ShellWrapper(snakemake).run(
    r"""
# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

if [[ "{args[filter_mode]}" != whole_genome ]]; then
    bedtools intersect -u -header -wa -a {snakemake.input.vcf} -b {path_bed} \
    | bcftools norm --remove-duplicates \
    | bcftools sort -o {snakemake.output.vcf} -O z
    tabix -f {snakemake.output.vcf}
else  # else, "all"
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
fi
"""
)
