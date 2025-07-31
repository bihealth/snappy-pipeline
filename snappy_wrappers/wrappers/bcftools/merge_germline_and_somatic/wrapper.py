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

The wrapper finds the loci of germline & somatic variants for one donor
and outputs them as a bed file, to use for CNV B-allele fraction computations.

Any germline variant overlapping with a somatic variant should be rejected,
but all somatic variants should be retained.
The user may want to select somatic variants after filtration.

``bcftools isec`` produces 4 files the temp directory. 3 of them are vcf files
containing somatic variants not in germline, germline variants not in somatic,
and common variants.
The fourth files is a tab-delimited file with 5 columns: the chromosome & position in
columns 1 & 2, the reference & alternate alleles in columns 3 & 4, and a flag in 
column 5, with value 10 when the locus in is germline only, 01 if the locus is somatic
only & 11 if the locus is in both.
The ``awk`` command transforms this text file into a bed file, rejecting overlapping loci
and indels (only SNVs are kept).
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

cmd = r"""
bcftools isec \
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
