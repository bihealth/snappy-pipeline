# -*- coding: utf-8 -*-
"""
Wrapper for the section of ``ascat.prepareHTS`` where the allele are counted.
Following usage example for version 3.2.0 (commit 3197356)

Mimicks C implementation of alleleCounter (https://github.com/cancerit/alleleCount)
which is required by ascat (https://github.com/VanLoo-lab/ascat), but not in bioconda nor conda-forge.

The replacement isn't perfect, and we are aware of two types of differences:
1. The loci not covered by any read at all are absent from the replacement script
2. Some variants have lower coverage in the original code. After extensive checks with samtools,
   we haven't been able to understand why this coverage should be lower. The output from the script
   seems in line with samtools mpileup.

Mandatory snakemake.input: bam, loci, reference
Optional snakemake.input:

Mandatory snakemake.params.args: max-depth, min-BQ, min-MQ, seed, skip-any-set,
    skip-any-unset
Optional snakemake.params.args:

Mandatory snakemake.output: freq, vcf
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
bcftools mpileup \
    --no-BAQ \
    --max-depth {args[max-depth]} \
    --skip-indels \
    --min-MQ {args[min-MQ]} \
    --min-BQ {args[min-BQ]} \
    --skip-any-set {args[skip-any-set]} --skip-any-unset {args[skip-any-unset]} \
    --fasta-ref "{snakemake.input.reference}" \
    --regions-file "{snakemake.input.loci}" \
    --seed {args[seed]} \
    "{snakemake.input.bam}" \
| bcftools annotate \
    --remove 'INFO,^FORMAT/AD' \
| bcftools view \
    --output-type z --output {snakemake.output.vcf} --write-index=tbi

# -------------------------------------------------------------------------------------------------
# Reformat vcf file to mimic files produced by `alleleCounter`:
# Tab-separated file with one row per SNP. The columns are:
# 1. Chromosome name
# 2. Position
# 3. Number of reads supporting an "A" at the position
# 4. Number of reads supporting an "C" at the position
# 5. Number of reads supporting an "G" at the position
# 6. Number of reads supporting an "T" at the position
# 7. Total number of reads at the position
#
bcftools query --format "%CHROM\t%POS\t%REF\t%ALT\t[%AD]" {snakemake.output.vcf} \
| awk -F '\t' '
BEGIN {{
    print "#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth"
}}
{{
    if ($3!~/^[ACGT]$/ || $4!~/^([ACGT],)*([ACGT]|<\*>)$/ || $3==$4) {{
        printf "Error in %s\n", $0
        exit 1
    }}

    a=0; c=0; g=0; t=0;

    patsplit($4, alts, /([ACGT]|<\*>)/);
    patsplit($5, counts, /[0-9]+/)

    if (length(alts)+1 != length(counts)) {{
        printf "Error in %s\n", $0
        exit 1
    }}

    if ($3=="A") {{ a += counts[1] }};
    if ($3=="C") {{ c += counts[1] }};
    if ($3=="G") {{ g += counts[1] }};
    if ($3=="T") {{ t += counts[1] }};
    
    s = counts[1];
    for (i=1; i<=length(alts); ++i) {{
        if (alts[i]=="A") {{ a += counts[i+1] }};
        if (alts[i]=="C") {{ c += counts[i+1] }};
        if (alts[i]=="G") {{ g += counts[i+1] }};
        if (alts[i]=="T") {{ t += counts[i+1] }};
        s += counts[i+1];
    }}
    
    printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, $2, a, c, g, t, s)
}}' \
> {snakemake.output.freq}
""".format(
    snakemake=snakemake,
    args=args,
)

ShellWrapper(snakemake).run(cmd)