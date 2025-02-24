# -*- coding: utf-8 -*-
"""
Wrapper for preparing allele counts for ASCAT
Following usage example for version 3.2.0 (commit 3197356)

`alleleCounter <https://github.com/cancerit/alleleCount>`_ is not available from conda-forge nor bioconda,
so the wrapper uses a like-for-like replacement shell script based on bcftools (``../alleleCounter/alleleCounter.sh``)
Note that to make ``bcftools`` happy, the genome sequence is added to the list of arguments,
which is not the case in the original script.

The ascat R script just dumps temporary files created by ``alleleCounter`` in the current directory.
To avoid the clutter, the ascat command is executed in a temp directory. Because of that, all paths to 
input & output files must be made absolute.

Mandatory snakemake.input: alleles, loci, normal, tumor
Optional snakemake.input: bed_file, probloci_file

Mandatory snakemake.params.args: chrom_names, gender, genomeVersion, minCounts, min_base_qual,
    min_map_qual, normalname, seed, tumorname
Optional snakemake.params.args: skip-any-set, skip-any-unset

Mandatory snakemake.output: normal_baf, normal_logr, tumor_baf, tumor_logr
Optional snakemake.output: 
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

PATTERN = re.compile(r"^(.*/)?(.+)_chr([0-9]+|[XY])\.txt$")

def find_prefix(snakemake_input):
    prefix = None
    for fn in map(str, snakemake_input):
        m = PATTERN.match(fn)
        if m:
            prefix = m.group(2) + "_chr"
            prefix = os.path.join(
                os.path.dirname(os.path.abspath(os.path.realpath(fn))),
                prefix
            )
            break
    assert prefix is not None, f"Can't find allele counts from {snakemake_input}"
    return prefix

allele_prefix = find_prefix(snakemake.input.alleles)

allele_frequencies_path = os.path.dirname(os.path.abspath(os.path.realpath(str(snakemake.input.tumor_alleleFrequencies[0]))))

cmd=r"""
library(ASCAT)

setwd("{allele_frequencies_path}")

# Note: normalseqfile not set to NA to avoid tumor only mode
ascat.prepareHTS(
    tumourseqfile=NA,
    normalseqfile="NA",
    skip_allele_counting_tumour=TRUE,
    skip_allele_counting_normal=TRUE,
    tumourname="{args[tumorname]}",
    normalname="{args[normalname]}",
    allelecounter_exe=NA,
    alleles.prefix="{allele_prefix}",
    gender="{args[gender]}",
    genomeVersion="{args[genomeVersion]}",
    nthreads={snakemake.resources._cores},
    tumourLogR_file="{tumourLogR_file}",
    tumourBAF_file="{tumourBAF_file}",
    normalLogR_file="{normalLogR_file}",
    normalBAF_file="{normalBAF_file}",
    minCounts={args[minCounts]},
    BED_file={bed_file},
    probloci_file={probloci_file},
    chrom_names={chrom_names},
    seed={args[seed]}
)
""".format(
    snakemake=snakemake,
    args=args,
    allele_frequencies_path=allele_frequencies_path,
    allele_prefix=allele_prefix,
    tumourLogR_file=os.path.abspath(os.path.realpath(str(snakemake.output.tumor_logr))),
    tumourBAF_file=os.path.abspath(os.path.realpath(str(snakemake.output.tumor_baf))),
    normalLogR_file=os.path.abspath(os.path.realpath(str(snakemake.output.normal_logr))),
    normalBAF_file=os.path.abspath(os.path.realpath(str(snakemake.output.normal_baf))),
    bed_file='"{}"'.format(os.path.abspath(os.path.realpath(str(snakemake.input.baits)))) if getattr(snakemake.input, "baits", None) is not None else "NA",
    probloci_file='"{}'.format(os.path.abspath(os.path.realpath(str(snakemake.input.probloci_file)))) if getattr(snakemake.input, "probloci_file", None) is not None else "NA",
    chrom_names='c("{}")'.format('", "'.join(args["chrom_names"])),
)

RWrapper(snakemake).run(cmd)