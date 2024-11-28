# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py reference"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = snakemake.params.get("args", {})

target = f"--target {args['target']}" if "target" in args else ""
antitarget = f"--antitarget {args['antitarget']}" if "antitarget" in args else ""

cmd = r"""
cnvkit.py reference \
    -o {snakemake.output.reference} \
    --fasta {args[reference]} \
    {cluster} {min_cluster_size} \
    {sample_sex} {male_reference} {diploid_parx_genome} \
    {no_gc} {no_edge} {no_rmask} \
    {target} {antitarget} {normals}
""".format(
    snakemake=snakemake,
    args=args,
    target=target,
    antitarget=antitarget,
    normals=" ".join(args["normals"]) if len(args.get("normals", [])) > 0 else "",
    cluster="--cluster" if args.get("cluster", False) else "",
    male_reference="--male-reference" if args.get("male-reference", False) else "",
    no_gc="--no-gc" if args.get("no-gc", False) else "",
    no_edge="--no-edge" if args.get("no-edge", False) else "",
    no_rmask="--no-rmask" if args.get("no-rmask", False) else "",
    min_cluster_size=f"--min-cluster-size {args['min-cluster-size']}" if args.get("min-cluster-size", None) is not None else "",
    sample_sex=f"--sample-sex {args['sample-sex']}" if args.get("sample-sex", None) is not None else "",
    diploid_parx_genome=f"--diploid-parx-genome {args['diploid-parx-genome']}" if args.get("diploid-parx-genome", None) is not None else ""
)

CnvkitWrapper(snakemake, cmd).run()
