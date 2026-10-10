# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing PureCN coverage"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

container = snakemake.input.container
intervals = snakemake.input.intervals

# Prepare files and directories that must be accessible by the container
files_to_bind = {
    "bam": snakemake.input.bam,
}

# Replace with full absolute paths
files_to_bind = {k: os.path.realpath(v) for k, v in files_to_bind.items()}
# Directories that mut be bound
dirs_to_bind = {k: os.path.dirname(v) for k, v in files_to_bind.items()}
# List of unique directories to bind: on cluster: <directory> -> from container: /bindings/d<i>)
bound_dirs = {e[1]: e[0] for e in enumerate(list(set(dirs_to_bind.values())))}
# Binding command
bindings = " ".join(["-B {}:/bindings/d{}:ro".format(k, v) for k, v in bound_dirs.items()])
# Path to files from the container
bound_files = {
    k: "/bindings/d{}/{}".format(bound_dirs[dirs_to_bind[k]], os.path.basename(v))
    for k, v in files_to_bind.items()
}

ShellWrapper(snakemake).run(
    r"""
# Create coverage
cmd="/usr/local/bin/Rscript /opt/PureCN/Coverage.R --force \
    --seed {config[seed]} \
    --out-dir $(dirname {snakemake.output.coverage}) \
    --bam {bound_files[bam]} \
    --intervals {intervals}
"
mkdir -p $(dirname {snakemake.output.coverage})
apptainer exec --home $PWD {bindings} {container} $cmd

# Rename coverage file name
d=$(dirname {snakemake.output.coverage})
bam_name=$(basename {snakemake.input.bam} .bam)
fn="$d/${{bam_name}}_coverage_loess.txt.gz"

test -e $fn
mv $fn {snakemake.output.coverage}
"""
)
