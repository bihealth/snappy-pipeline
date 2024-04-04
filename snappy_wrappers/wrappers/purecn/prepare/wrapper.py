# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing exome kit intervals for PureCN"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

genome = snakemake.config["static_data_config"]["reference"]["path"]

# Prepare files and directories that must be accessible by the container
bound_files = {
    "genome": os.path.normpath(genome),
    "path_bait_regions": os.path.normpath(config["path_bait_regions"]),
    "mappability": os.path.normpath(config["mappability"])
    if "mappability" in config and config["mappability"]
    else "",
    "reptiming": os.path.normpath(config["reptiming"])
    if "reptiming" in config and config["reptiming"]
    else "",
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

md5() {{
    d=$(dirname $1)
    f=$(basename $1)
    pushd $d 1> /dev/null 2>&1
    md5sum $f > $f.md5
    popd 1> /dev/null 2>&1
}}

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

md5 {snakemake.output.intervals}
md5 {snakemake.output.optimized}
md5 {snakemake.output.optimized}.tbi
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
