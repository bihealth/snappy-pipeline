# -*- coding: utf-8 -*-
"""
Wrapper for preparing allele counts for ASCAT

Mandatory snakemake.input: alleles, tumorAlleleCounts
Optional snakemake.input: bed_file, normalAlleleCounts, probloci_file

Mandatory snakemake.params.args: chrom_names, gender, genomeVersion, min_base_qual, minCounts,
    min_map_qual, normalname, seed, tumorname
Optional snakemake.params.args:

Mandatory snakemake.output: tumor_baf, tumor_logr
Optional snakemake.output: normal_baf, normal_logr
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import re
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.snappy_wrapper import RWrapper

args = getattr(snakemake.params, "args", {})

has_normal = getattr(snakemake.output, "normal_logr", None) is not None and getattr(snakemake.output, "normal_baf", None) is not None

PATTERN = re.compile(r"^(.+)_chr([0-9]+|[XY])\.txt$")

def find_prefix(snakemake_input):
    prefix = None
    for fn in map(str, snakemake_input):
        m = PATTERN.match(fn)
        if m:
            prefix = m.groups()[0] + "_chr"
            break
    assert prefix is not None, f"Can't find allele counts from {snakemake_input}"
    return prefix

allele_prefix = find_prefix(snakemake.input.alleles)
tumorAlleleCounts_prefix = find_prefix(snakemake.input.tumorAlleleCounts)
if has_normal:
    normalAlleleCounts_prefix = '"' + find_prefix(snakemake.input.normalAlleleCounts) + '"'
else:
    normalAlleleCounts_prefix = "NA"

cmd=r"""
library(ASCAT)

# Takes the part of ``ascat.prepareHTS`` which doesn't invoke ``allelecounter``.
# That is: ``ascat.getBAFsAndLogRs`` & ``ascat.synchroniseFiles``

ascat.getBAFsAndLogRs(
    samplename="{args[tumorname]}",
    tumourAlleleCountsFile.prefix="{tumorAlleleCounts_prefix}",
    normalAlleleCountsFile.prefix={normalAlleleCounts_prefix},
    tumourLogR_file="{snakemake.output.tumor_logr}",
    tumourBAF_file="{snakemake.output.tumor_baf}",
    normalLogR_file={normal_logr},
    normalBAF_file={normal_baf},
    alleles.prefix="{allele_prefix}",
    gender="{args[gender]}",
    genomeVersion="{args[genomeVersion]}",
    chrom_names={chrom_names},
    minCounts={args[minCounts]},
    BED_file={bed_file},
    probloci_file={probloci_file},
    tumour_only_mode={tumor_only},
    loci_binsize=1,
    seed={args[seed]}
)

ascat.synchroniseFiles(
    samplename="{args[tumorname]}",
    tumourLogR_file="{snakemake.output.tumor_logr}",
    tumourBAF_file="{snakemake.output.tumor_baf}",
    normalLogR_file={normal_logr},
    normalBAF_file={normal_baf}
)
""".format(
    snakemake=snakemake,
    args=args,
    tumorAlleleCounts_prefix=tumorAlleleCounts_prefix,
    normalAlleleCounts_prefix=normalAlleleCounts_prefix,
    allele_prefix=allele_prefix,
    normal_logr=f'"{snakemake.output.normal_logr}"' if has_normal else "NA",
    normal_baf=f'"{snakemake.output.normal_baf}"' if has_normal else "NA",
    bed_file=f'"{snakemake.input.bed_file}"' if getattr(snakemake.input, "bed_file", None) is not None else "NA",
    probloci_file=f'"{snakemake.input.probloci_file}"' if getattr(snakemake.input, "probloci_file", None) is not None else "NA",
    chrom_names='c("{}")'.format('", "'.join(args["chrom_names"])),
    tumor_only="FALSE" if has_normal else "TRUE",
)

RWrapper(snakemake).run(cmd)