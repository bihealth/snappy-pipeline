# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing PureCN coverage
"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

if "container" in snakemake.input.keys() and snakemake.input.container:
    container = snakemake.input.container
elif "path_container" in config.keys() and config["path_container"]:
    container = config["path_container"]
else:
    raise Exception("No path to PureCN container")

if "intervals" in snakemake.input.keys() and snakemake.input.intervals:
    intervals = snakemake.input.intervals
elif "path_intervals" in config.keys() and config["path_intervals"]:
    intervals = config["path_intervals"]
else:
    raise Exception("No path to PureCN intervals")

# Prepare files and directories that must be accessible by the container
files_to_bind = {
    "bam": snakemake.input.bam,
}
if "intervals" not in snakemake.input.keys():
    files_to_bind["intervals"] = intervals

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

if "intervals" in bound_files.keys():
    intervals = bound_files["intervals"]

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

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
mapper="{snakemake.wildcards[mapper]}"
libname="{snakemake.wildcards[library_name]}"
fn="$d/$mapper.${{libname}}_coverage_loess.txt.gz"

test -e $fn
mv $fn {snakemake.output.coverage}

pushd $(dirname {snakemake.output.coverage})
md5sum $(basename {snakemake.output.coverage}) > $(basename {snakemake.output.coverage}).md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
