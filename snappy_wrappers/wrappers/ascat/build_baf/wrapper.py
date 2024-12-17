# -*- coding: utf-8 -*-
"""
Wrapper for building BAF files for ASCAT.

It is a replacement for ``alleleCounter`` (which is not part of the bioconda package)

Mandatory snakemake.input: bam, locii, reference
Optional snakemake.input: baits

Mandatory snakemake.params.args: chrom_names, exclude_flags, gender, genomeVersion, include_flags,
    max_coverage, min_base_qual, minCounts, min_map_qual
Optional snakemake.params.args:

Mandatory snakemake.output: txt, vcf
Optional snakemake.output:
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

if getattr(snakemake.input, "baits", None) is not None:
    prefix=r"""
        if [[ {snakemake.input.baits} =~ \.gz$ ]]
        then
            zcat {snakemake.input.baits} | bgzip > $TMPDIR/baits.bed.gz
        else
            bgzip {snakemake.input.baits} > $TMPDIR/baits.bed.gz
        fi
        tabix $TMPDIR/baits.bed.gz

    """
    baits = "-R $TMPDIR/baits.bed.gz"
else:
    prefix = ""
    baits = ""

cmd=r"""
# Perform pileups at the spot positions.
# ASCAT default values:
# min_mapq: 35, min_baq: 20, exclude_flags: 3582, include_flags: 3, min_depth: 10
#
bcftools mpileup \
    -B \
    -d {args[max_coverage]} \
    -I \
    -q {args[min_map_qual]} \
    -Q {args[min_base_qual]} \
    --ns {args[exclude_flags]} --nu {args[include_flags]} \
    -f {snakemake.input.reference} \
    -R {snakemake.input.locii} \
    {snakemake.input.bam} \
| bcftools annotate \
    -x 'INFO,^FORMAT/AD' \
| bcftools view \
    -A \
    -i 'FORMAT/AD[0:0] + FORMAT/AD[0:1] > {args[minCounts]}' \
    {baits} \
    -O z -o {snakemake.output.vcf} -W

# -------------------------------------------------------------------------------------------------
# Reformat vcf file to mimic files produced by `alleleCounter`
#
zgrep -v '^#' {snakemake.output.vcf} \
| cut -f 1-2,4-5,10 \
| awk -F '\t' '
BEGIN {{
    print "#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth"
}}
{{
    if ($3!~/^[ACGT]$/ || $4!~/^([ACGT](,[ACGT])*|<\*>)$/ || $3==$4) {{
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
> {snakemake.output.txt}
""".format(snakemake=snakemake, args=args, baits=baits)

ShellWrapper(snakemake).run(prefix + cmd)