# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py scatter"""

import csv
import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

# Fix chromosome name prefix
if args.get("chromosome", None) is not None:
    chromosome = args["chromosome"]
    if chromosome.startswith("chr"):
        ucsc = chromosome
        ensembl = chromosome[3:]
        if ensembl == "M":
            ensembl = "MT"
    else:
        ucsc = f"chr{chromosome}"
        ensembl = chromosome
        if ucsc == "chrMT":
            ucsc = "chrM"

    with open(snakemake.input.segments, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for record in reader:
            if ucsc == record["chromosome"]:
                args["chromosome"] = ucsc
                break
            if ensembl == record["chromosome"]:
                args["chromosome"] = ensembl
                break

if snakemake.input.get("variants", None) is not None:
    variants = r"""
        ---vcf {snakemake.input.variants} \
        --sample-id {args[sample-id]} --normal-id {args[normal-id]} \
        --min-variant-depth {args[min-variant-depth]} {zygocity_freq}
    """.format(
        snakemake=snakemake,
        args=args,
        zygocity_freq=f"--zygocity-freq {args['zygocity-freq']}" if args.get("zygocity-freq", None) is not None else "",
    )
else:
    variants = ""

cmd = r"""
cnvkit.py scatter \
    -o {snakemake.output.plot} \
    --segment {snakemake.input.segments} \
    {chromosome} {gene} {range_list} \
    --width {args[width]} \
    --antitarget-marker {args[antitarget-marker]} --segment-color {args[segment-color]} \
    {by_bin} {trend} --title "{args[title]}" \
    {y_min} {y_max} {fig_size} \
    {variants} \
    {snakemake.input.ratios}
""".format(
    snakemake=snakemake,
    args=args,
    variants=variants,
    chromosome=f"--chromosome {args['chromosome']}" if args.get("chromosome", None) is not None else "",
    gene=f"--gene {args['gene']}" if args.get("gene", None) is not None else "",
    range_list=f"--range-list {snakemake.input.range_list}" if snakemake.input.get("range_list", None) is not None else "",
    by_bin="--by-bin" if args.get("by-bin", False) else "",
    trend="--trend" if args.get("trend", False) else "",
    y_min=f"--y-min {args['y-min']}" if args.get("y-min", None) is not None else "",
    y_max=f"--y-max {args['y-max']}" if args.get("y-max", None) is not None else "",
    fig_size="--fig-size {}".format(" ".join(map(str, args['fig-size']))) if args.get("fig-size", None) is not None else "",
)

CnvkitWrapper(snakemake, cmd).run()
