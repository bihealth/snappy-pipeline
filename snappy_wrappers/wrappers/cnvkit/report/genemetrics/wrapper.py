# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py genemetrics"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

cmd = r"""
cnvkit.py genemetrics \
    -o {snakemake.output.report} \
    --segment {snakemake.input.segments} \
    --threshold {args[threshold]} --min-probes {args[min-probes]} \
    {drop_low_coverage} {male_reference} {sample_sex} {diploid_parx_genome} \
    {stats} \
    {snakemake.input.ratios}
""".format(
    snakemake=snakemake,
    args=args,
    drop_low_coverage="--drop-low-coverage" if args.get("drop-low-coverage", False) else "",
    male_reference="--male-reference" if args.get("male-reference", False) else "",
    sample_sex=f"--sample-sex {args['sample-sex']}" if args.get("sample-sex", None) is not None else "",
    diploid_parx_genome=f"--diploid-parx-genome {args['diploid-parx-genome']}" if args.get("diploid-parx-genome", None) is not None else "",
    stats=" ".join([f"--{stat}" for stat in args["stats"]]),
)

CnvkitWrapper(snakemake, cmd).run()
