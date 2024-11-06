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

cmd = r"""
cnvkit.py reference \
    -o {snakemake.output.panel} \
    --fasta {args[reference]} \
    {cluster} {min_cluster_size} \
    {sample_sex} {male_reference} {diploid_parx_genome} \
    {no_gc} {no_edge} {no_rmask} \
    {target} {antitarget} {normals}
""".format(
    snakemake=snakemake,
    args=args,
    cluster="--cluster" if "cluster" in args else "",
    min_cluster_size=f"--min-cluster-size {args['min_cluster_size']}" if "cluster" in args and "min_cluster_size" in args else "",
    no_gc="--no-gc" if "no_gc" in args else "",
    no_edge="--no-edge" if "no_edge" in args else "",
    no_rmask="--no-rmask" if "no_rmask" in args else "",
    sample_sex=f"--sample-sex {args['sample_sex']}" if "sample_sex" in args else "",
    male_reference="--male-reference" if "male_reference" in args else "",
    diploid_parx_genome=f"--diploid_parx_genome {args['diploid_parx_genome']}" if "diploid_parx_genome" in args else "",
    target=f"--target {snakemake.input.target}" if snakemake.input.get("target", None) else "",
    antitarget=f"--antitarget {snakemake.input.antitarget}" if  snakemake.input.get("antitarget", None) else "",
    normals=" ".join(snakemake.input.normals) if snakemake.input.get("normals", None) else "",
)

CnvkitWrapper(snakemake, cmd).run()
