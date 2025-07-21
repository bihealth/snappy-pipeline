# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for applying the filter list."""

import gzip
import re
import sys
import textwrap

from snakemake import shell

args = getattr(snakemake.params, "args", {})
config = args.get("config", {})
print("DEBUG- args = {}".format(args), file=sys.stderr)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

title = re.compile(
    "^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER(\tINFO)?(\tFORMAT)?\t([^\t]+)\t([^\t]+)"
)
iTumor = -1
with gzip.open(snakemake.input.vcf, mode="rt") as f:
    for line in f:
        m = title.search(line.rstrip())
        if m:
            if (m.group(3) == args["normal_sample"]) and (m.group(4) == args["tumor_sample"]):
                iTumor = 1
                break
            if (m.group(4) == args["normal_sample"]) and (m.group(3) == args["tumor_sample"]):
                iTumor = 0
                break
if iTumor < 0:
    print(
        "Can't find normal sample {} or tumor sample {} in vcf header".format(
            args["normal_sample"], args["tumor_sample"]
        ),
        file=sys.stderr,
    )
    sys.exit(-1)

cmd = ["zcat " + snakemake.input.vcf]

if "dkfz" in args["filter_set"]:
    cmd.append("bcftools view -f .,PASS")
if "ebfilter" in args["filter_set"]:
    cmd.append(
        "bcftools view -e 'EB < {threshold}'".format(
            threshold=config["dkfz_and_ebfilter"].get("ebfilter_threshold", 3)
        )
    )
if "oxog" in args["filter_set"]:
    if args["var_caller"] != "scalpel":
        min_vaf = config["dkfz_and_ebfilter_and_oxog"].get("vaf_threshold", 0.08)
        min_cov = config["dkfz_and_ebfilter_and_oxog"].get("coverage_threshold", 5)
        if args["var_caller"] == "mutect2":
            allele_freq_str = "AF"
        else:
            allele_freq_str = "FA"
        oxo_filter = '(FILTER != "PASS" || FORMAT/{}[{}:0]<={} || FORMAT/AD[{}:1]<{})'.format(
            allele_freq_str, iTumor, min_vaf, iTumor, min_cov
        )
        cmd.append(
            'bcftools filter -e \'{oxo_filter} && REF = "G" && ALT ~ "T"\''.format(
                oxo_filter=oxo_filter
            )
        )
        cmd.append(
            'bcftools filter -e \'{oxo_filter} && REF = "C" && ALT ~ "A"\''.format(
                oxo_filter=oxo_filter
            )
        )

script = " | ".join(cmd) + " | bgzip > " + snakemake.output.vcf

shell(
    textwrap.dedent(
        r"""
    set -x

    {script}

    tabix {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi}).md5
        """
    )
)

# TODO
# Compute MD5 sums of logs.
# shell(
#    r"""
# pwd
# md5sum {snakemake.log} >{snakemake.log}.md5
# """
# )
