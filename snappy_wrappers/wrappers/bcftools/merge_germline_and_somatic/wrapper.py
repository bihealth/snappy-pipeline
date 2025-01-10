# -*- coding: utf-8 -*-
"""
Wrapper to merge germline & somatic SNVs locii into one location file (bed)

It is meant to be used in conjunction with other bcftools commands, such as mpileup & call

Mandatory snakemake.input: germline, somatic
Optional snakemake.input: regions

Mandatory snakemake.params.args:
Optional snakemake.params.args:

Mandatory snakemake.output: bed
Optional snakemake.output:
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

cmd = r"""
bcftools iseq \
    --collapse all --regions-overlap 2 \
    {regions} \
    --prefix $TMPDIR/ \
    {snakemake.input.germline} {snakemake.input.somatic}

awk \
    -F '\t' \
    '($5 != "11") && (length($3) == 1) && (length($3) == length($4)) {{printf "%s\t%d\t%d\n", $1, $2-1, $2}}' \
    $TMPDIR/sites.txt \
| bgzip \
> {snakemake.output.bed}

tabix {snakemake.output.bed}
""".format(
    snakemake=snakemake,
    regions=f"--regions-file {snakemake.input.regions}" if getattr(snakemake.input, regions, None) is not None else "",
)

ShellWrapper(snakemake).run(cmd)
