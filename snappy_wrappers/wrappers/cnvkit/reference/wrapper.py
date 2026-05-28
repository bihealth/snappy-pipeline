# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py reference"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

# NOTE: snakemake.input.target and snakemake.input.antitarget contain
#       the output of target & antitarget substeps when there is no bam files
#       the bam files lists when the list of normals is not empty

cluster = (
    " --cluster --min-cluster-size {}".format(args["min_cluster_size"])
    if args["min_cluster_size"] > 0
    else ""
)
gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""
no_gc = " --no-gc" if not args["gc_correction"] else ""
no_edge = " --no-edge" if not args["edge_correction"] else ""
no_rmask = " --no-rmask" if not args["rmask_correction"] else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

if [[ "{args[flat]}" = "True" ]]
then
    cnvkit.py reference \
        --output {snakemake.output.panel} \
        --fasta {snakemake.input.reference} \
        {cluster} {gender} {male} {no_gc} {no_edge} {no_rmask} \
        --targets {snakemake.input.target} --antitargets {snakemake.input.antitarget}
else
    cnvkit.py reference \
        --output {snakemake.output.panel} \
        --fasta {snakemake.input.reference} \
        {cluster} {gender} {male} {no_gc} {no_edge} {no_rmask} \
        {snakemake.input.target} {snakemake.input.antitarget}
fi

if [[ -n "{snakemake.input.logs}" ]]
then
    tar -zcvf {snakemake.output.log} {snakemake.input.logs}
else
    touch {snakemake.output.log}
fi
"""
)
