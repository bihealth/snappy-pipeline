# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py reference"""

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

cmd = r"""
cnvkit.py reference \
    -o {snakemake.output.reference} \
    --fasta {snakemake.params.reference} \
    {cluster} {min_cluster_size} \
    {sample_sex} {male_reference} {diploid_parx_genome} \
    {no_gc} {no_edge} {no_rmask} \
    {target} {antitarget} {normals}
""".format(
    snakemake=snakemake,
    cluster="--cluster" if snakemake.params.cluster else "",
    min_cluster_size=f"--min-cluster-size {snakemake.params.min_cluster_size}" if snakemake.params.cluster and "min_cluster_size" in snakemake.params else "",
    no_gc="--no-gc" if snakemake.params.no_gc else "",
    no_edge="--no-edge" if snakemake.params.no_edge else "",
    no_rmask="--no-rmask" if snakemake.params.no_rmask else "",
    sample_sex=f"--sample-sex {snakemake.params.sample_sex}" if "sample_sex" in snakemake.params else "",
    male_reference="--male-reference" if snakemake.params.male_reference else "",
    diploid_parx_genome=f"--diploid_parx_genome {snakemake.params.diploid_parx_genome}" if "diploid_parx_genome" in snakemake.params else "",
    target=f"--target {snakemake.input.target}" if "target" in snakemake.input else "",
    antitarget=f"--antitarget {snakemake.input.antitarget}" if "antitarget" in snakemake.input else "",
    normals=" ".join(snakemake.input.normals) if "normals" in snakemake.input else "",
)

CnvkitWrapper(snakemake, cmd).run()
