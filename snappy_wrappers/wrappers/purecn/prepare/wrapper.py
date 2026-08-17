# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing exome kit intervals for PureCN"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["purecn"]

genome = snakemake.input.reference

# Prepare files and directories that must be accessible by the container
bound_files = {
    "genome": os.path.normpath(genome),
    "path_bait_regions": os.path.normpath(config["path_bait_regions"]),
    "mappability": (
        os.path.normpath(config["mappability"])
        if "mappability" in config and config["mappability"]
        else ""
    ),
    "reptiming": (
        os.path.normpath(config["reptiming"])
        if "reptiming" in config and config["reptiming"]
        else ""
    ),
}

keys = list(bound_files.keys())
bindings = []
for i in range(len(keys)):
    k = keys[i]
    if bound_files[k]:
        # Binding directory to /bindings/d<i>
        bindings.append(" -B {}:/bindings/d{}:ro".format(os.path.dirname(bound_files[k]), i))
        bound_files[k] = "/bindings/d{}/{}".format(i, os.path.basename(bound_files[k]))
bindings = " ".join(bindings)

ShellWrapper(snakemake).run(
    r"""
# Create panel
uncompressed=$(echo "{snakemake.output.optimized}" | sed -e "s/\.gz$//")

cmd="/usr/local/bin/Rscript /opt/PureCN/IntervalFile.R \
    --out-file {snakemake.output.intervals} \
    --export $uncompressed \
    --in-file {bound_files[path_bait_regions]} \
    --fasta {bound_files[genome]} --genome {config[genome_name]} \
    $(if [[ -n "{bound_files[mappability]}" ]]; then \
        echo "--mappability {bound_files[mappability]}"
    fi) \
    $(if [[ -n "{bound_files[reptiming]}" ]]; then \
        echo "--reptiming {bound_files[reptiming]}"
    fi)
"

mkdir -p $(dirname {snakemake.output.intervals})
mkdir -p $(dirname $uncompressed)

apptainer exec --home $PWD {bindings} {snakemake.input.container} $cmd

bgzip $uncompressed
tabix {snakemake.output.optimized}
"""
)
